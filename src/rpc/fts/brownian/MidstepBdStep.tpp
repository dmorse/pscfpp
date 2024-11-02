#ifndef RPC_MIDSTEP_BD_STEP_TPP
#define RPC_MIDSTEP_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MidstepBdStep.h"

#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/System.h>
#include <pscf/math/IntVec.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   MidstepBdStep<D>::MidstepBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      wh_(),
      wf_(),
      eta_(),
      dwc_(),
      dwp_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   MidstepBdStep<D>::~MidstepBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void MidstepBdStep<D>::readParameters(std::istream &in)
   {
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      wf_.allocate(nMonomer);
      wh_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         wf_[i].allocate(meshDimensions);
         wh_[i].allocate(meshDimensions);
      }
      eta_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         eta_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);
      dwp_.allocate(meshDimensions);

   }

   template <int D>
   void MidstepBdStep<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(wf_.capacity() == nMonomer);
      UTIL_CHECK(wh_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(wf_[i].capacity() == meshSize);
         UTIL_CHECK(wh_[i].capacity() == meshSize);
      }
      UTIL_CHECK(eta_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(eta_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      UTIL_CHECK(dwp_.capacity() == meshSize);
   }

   template <int D>
   bool MidstepBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Save current state
      simulator().saveState();

      // Copy current W fields from parent system
      for (i = 0; i < nMonomer; ++i) {
         wh_[i] = system().w().rgrid(i);
         wf_[i] = wh_[i];
      }

      // Store initial value of pressure field
      dwp_ = simulator().wc(nMonomer-1);

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);

      // Construct all random displacements
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> & eta = eta_[j];
         for (k = 0; k < meshSize; ++k) {
            eta[k] = b*random().gaussian();
         }
      }

      // Mid-step Predictor
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         RField<D> const & eta = eta_[j];
         for (k = 0; k < meshSize; ++k) {
            dwc_[k] = a*dc[k] + eta[k];
         }

         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wn = wh_[i];
            evec = 0.5*simulator().chiEvecs(j,i);
            for (k = 0; k < meshSize; ++k) {
               wn[k] += evec*dwc_[k];
            }
         }

      }

      // Set modified fields at mid-point
      system().setWRGrid(wh_);

      // Enforce incompressibility (also solves MDE repeatedly)
      bool isConverged = false;
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else{
         UTIL_CHECK(system().hasCFields());

         // Compute components and derivatives at mid-point
         simulator().clearData();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();

         // Store change in pressure field
         RField<D> const & wp = simulator().wc(nMonomer-1);
         for (k = 0; k < meshSize; ++k) {
            dwp_[k] = wp[k] - dwp_[k];
         }

         // Full step (corrector)
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D> const & dc = simulator().dc(j);
            RField<D> const & eta = eta_[j];
            for (k = 0; k < meshSize; ++k) {
               dwc_[k] = a*dc[k] + eta[k];
            }
            for (i = 0; i < nMonomer; ++i) {
               RField<D> & wn = wf_[i];
               evec = simulator().chiEvecs(j,i);
               for (k = 0; k < meshSize; ++k) {
                  wn[k] += evec*dwc_[k];
               }
            }
         }

         // Extrapolate pressure field
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wn = wf_[i];
            for (k = 0; k < meshSize; ++k) {
               wn[k] += 2.0*dwp_[k];
            }
         }

         // Set fields at final point
         system().setWRGrid(wf_);

         int compress2 = simulator().compressor().compress();
         if (compress2 != 0){
            simulator().restoreState();
         } else {
            isConverged = true;
            UTIL_CHECK(system().hasCFields());

            // Compute components and derivatives at final point
            simulator().clearState();
            simulator().clearData();
            simulator().computeWc();
            simulator().computeCc();
            simulator().computeDc();

         }
      }
      return isConverged;
   }

}
}
#endif
