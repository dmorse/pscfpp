#ifndef PSPG_EXPLICIT_BD_STEP_TPP
#define PSPG_EXPLICIT_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExplicitBdStep.h"

#include <pspg/simulate/bdstep/BdSimulator.h>
#include <pspg/compressor/Compressor.h>
#include <pspg/System.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   ExplicitBdStep<D>::ExplicitBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      w_(),
      dwc_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   ExplicitBdStep<D>::~ExplicitBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void ExplicitBdStep<D>::readParameters(std::istream &in)
   {
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);

   }

   template <int D>
   void ExplicitBdStep<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      gaussianField_.allocate(meshSize);
      dwd_.allocate(meshSize);
   }

   template <int D>
   void ExplicitBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      // Copy current W fields from parent system into w_
      DArray<RField<D>> const * Wr = &system().w().rgrid();
      for (i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (w_[i].cField(), (*Wr)[i].cField(), meshSize);
      }
      
      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      double a = -1.0*mobility_;
      double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);
      // Constants for normal distribution
      double stddev = 1.0;
      double mean = 0;
      
      // Modify local field copy wc_
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         // Generagte normal distributed random floating point numbers
         cudaRandom().normal(gaussianField_.cField(), meshSize, (cudaReal)stddev, (cudaReal)mean);
         // dwr
         scaleReal<<<nBlocks, nThreads>>>(gaussianField_.cField(), b, meshSize);
         // dwd
         assignReal<<<nBlocks, nThreads>>>(dwd_.cField(), dc.cField(), meshSize);
         scaleReal<<<nBlocks, nThreads>>>(dwd_.cField(), a, meshSize);
         // dwc
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>
            (dwd_.cField(), gaussianField_.cField(), dwc_.cField(), meshSize);
         
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            pointWiseAddScale<<<nBlocks, nThreads>>>
               (w.cField(), dwc_.cField(), evec, meshSize);
         }
      }

      // Set modified fields in parent system
      system().setWRGrid(w_);
      simulator().clearData();

      // Enforce incompressibility (also solves MDE repeatedly)
      system().compressor().compress();
      UTIL_CHECK(system().hasCFields());

      // Evaluate component properties in new state
      simulator().computeWc();
      simulator().computeCc();
      simulator().computeDc();
   }

}
}
#endif
