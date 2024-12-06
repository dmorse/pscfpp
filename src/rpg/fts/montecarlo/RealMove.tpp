#ifndef RPG_REAL_MOVE_TPP
#define RPG_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h" 
#include <pscf/math/IntVec.h>
#include <util/param/ParamComposite.h>
#include <rpg/System.h>

namespace Pscf {
namespace Rpg {

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
      const IntVec<D> dimensions = system().domain().mesh().dimensions();
      if (!isAllocated_){
         wFieldTmp_.allocate(nMonomer);
         randomField_.allocate(dimensions);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
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
      Array<RField<D>> const * currSys = &system().w().rgrid();

      // For multi-component copolymer
      for (int i = 0; i < nMonomer; i++){

         // Generate random numbers between 0.0 and 1.0 from uniform distribution
         cudaRandom().uniform(randomField_.cArray(), meshSize);

         // Generate random numbers between [-stepSize_,stepSize_]
         mcftsScale<<<nBlocks, nThreads>>>(randomField_.cArray(), stepSize_, meshSize);

         // Change the w field configuration
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>((*currSys)[i].cArray(), randomField_.cArray(),
                                                      wFieldTmp_[i].cArray(), meshSize);

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
   
   template<int D>
   void RealMove<D>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << "RealMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif
