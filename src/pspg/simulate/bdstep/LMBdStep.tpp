#ifndef PSPG_LM_BD_STEP_TPP
#define PSPG_LM_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LMBdStep.h"

#include <pspg/simulate/bdstep/BdSimulator.h>
#include <pspg/compressor/Compressor.h>
#include <pspg/System.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/CudaRandom.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   LMBdStep<D>::LMBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      w_(),
      etaA_(),
      etaB_(),
      dwc_(),
      etaNewPtr_(0),
      etaOldPtr_(0),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   LMBdStep<D>::~LMBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void LMBdStep<D>::readParameters(std::istream &in)
   {
      read(in, "mobility", mobility_);

      int nMonomer = system().mixture().nMonomer();
      int meshSize = system().domain().mesh().size();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();

      // Allocate memory for private containers
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      etaA_.allocate(nMonomer-1);
      etaB_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         etaA_[i].allocate(meshDimensions);
         etaB_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);
      gaussianField_.allocate(meshSize);
   }

   /*
   * Generate new random displacement values
   */
   template <int D>
   void LMBdStep<D>::generateEtaNew()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      // GPU resources  
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Prefactor b for displacements
      const double vSystem = system().domain().unitCell().volume();
      double b = sqrt(0.5*mobility_*double(meshSize)/vSystem);
      // Constants for normal distribution
      double stddev = 1.0;
      double mean = 0;

      int j;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D>& eta = etaNew(j);
         // Generagte normal distributed random floating point numbers
         cudaRandom().normal(gaussianField_.cField(), meshSize, (cudaReal)stddev, (cudaReal)mean);
         scaleReal<<<nBlocks, nThreads>>>(gaussianField_.cField(), b, meshSize);
         assignReal<<<nBlocks, nThreads>>>(eta.cField(), gaussianField_.cField(), meshSize);
      }
   }

   template <int D>
   void LMBdStep<D>::exchangeOldNew()
   {
      DArray< RField<D> >* temp;
      temp = etaOldPtr_;
      etaOldPtr_ = etaNewPtr_;
      etaNewPtr_ = temp;
   }

   /*
   * Initial setup before main simulation loop.
   */
   template <int D>
   void LMBdStep<D>::setup()
   {
      int nMonomer = system().mixture().nMonomer();
      int meshSize = system().domain().mesh().size();

      // Check array capacities
      UTIL_CHECK(etaA_.capacity() == nMonomer-1);
      UTIL_CHECK(etaB_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(etaA_[i].capacity() == meshSize);
         UTIL_CHECK(etaB_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);

      // Initialize pointers
      etaOldPtr_ = &etaA_;
      etaNewPtr_ = &etaB_;

      generateEtaNew();
      exchangeOldNew();
   }

   /*
   * One step of Leimkuhler-Matthews BD algorithm.
   */
   template <int D>
   void LMBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Copy current fields to w_
      for (i = 0; i < nMonomer; ++i) {
         w_[i] = system().w().rgrid(i);
      }

      // Generate new random displacement values
      generateEtaNew();

      // Take LM step:
      double a = -1.0*mobility_;
      double evec;
      // Loop over composition eigenvectors of projected chi matrix
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & etaN = etaNew(j);
         RField<D> const & etaO = etaOld(j);
         RField<D> const & dc = simulator().dc(j);
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>
            (etaN.cField(), etaO.cField(), dwc_.cField(), meshSize);
         pointWiseAddScale<<<nBlocks, nThreads>>>
            (dwc_.cField(), dc.cField(), a, meshSize);
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            pointWiseAddScale<<<nBlocks, nThreads>>>(w.cField(), dwc_.cField(), evec, meshSize);
         }
      }

      // Set modified fields at predicted state wp_
      system().setWRGrid(w_);
      system().compressor().compress();
      UTIL_CHECK(system().hasCFields());

      // Compute components and derivatives at wp_
      simulator().clearData();
      simulator().computeWc();
      simulator().computeCc();
      simulator().computeDc();

      // Exchange old and new random fields
      exchangeOldNew();
   }

}
}
#endif
