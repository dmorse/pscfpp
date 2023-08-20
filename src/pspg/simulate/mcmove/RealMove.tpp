#ifndef PSPG_REAL_MOVE_TPP
#define PSPG_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h" 
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <pspg/System.h>
#include <util/archives/Serializable_includes.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <sys/time.h>


namespace Pscf {
namespace Pspg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealMove<D>::RealMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      isAllocated_(false)
   { setClassName("RealMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   RealMove<D>::~RealMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void RealMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", stepSize_);
   }
   
   template <int D>
   void RealMove<D>::setup()
   {  
      McMove<D>::setup();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      if (!isAllocated_){
         wFieldTmp_.allocate(nMonomer);
         randomField_.allocate(meshSize);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(meshSize);
         }
         isAllocated_ = true;
      }
      // Create pseudo-random number generator on gpu
      CURAND_CALL(curandCreateGenerator(&gen_, CURAND_RNG_PSEUDO_DEFAULT));
      // Set seed
      unsigned long long seed;
      timeval time;
      gettimeofday(&time, NULL);
      seed = time.tv_sec + 1123*time.tv_usec;
      //CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen_, seed));
      curandStatus_t status;
      status = curandSetPseudoRandomGeneratorSeed(gen_, seed);
      if(status != CURAND_STATUS_SUCCESS){
         std::cout<<"Generator random number error "<<std::endl;
      }
      
   }
   
   
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void RealMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      Array<RDField<D>> const * currSys = &system().w().rgrid();
      
      // For multi-component copolymer
      for (int i = 0; i < nMonomer; i++){
         // Generate random numbers between 0.0 and 1.0 from uniform distribution
#ifdef SINGLE_PRECISION
         curandStatus_t gen_error = curandGenerateUniform(gen_, randomField_.cDField(), meshSize); 
#else
         curandStatus_t gen_error = curandGenerateUniformDouble(gen_, randomField_.cDField(), meshSize);
#endif
         // Generate random numbers between [-stepSize_,stepSize_]
         mcftsScale<<<nBlocks, nThreads>>>(randomField_.cDField(), stepSize_, meshSize);
         
         // Change the w field configuration
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>((*currSys)[i].cDField(), randomField_.cDField(), 
                                                      wFieldTmp_[i].cDField(), meshSize);
         
      }
      // set system r grid
      system().setWRGrid(wFieldTmp_);

   }


   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void RealMove<D>::output()
   {}

}
}
#endif
