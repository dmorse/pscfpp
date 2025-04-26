#ifndef RPC_PRED_CORR_BD_STEP_TPP
#define RPC_PRED_CORR_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"

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
   PredCorrBdStep<D>::PredCorrBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      wp_(),
      wf_(),
      dci_(),
      eta_(),
      dwc_(),
      dwp_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   PredCorrBdStep<D>::~PredCorrBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void PredCorrBdStep<D>::readParameters(std::istream &in)
   {
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      wp_.allocate(nMonomer);
      wf_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         wp_[i].allocate(meshDimensions);
         wf_[i].allocate(meshDimensions);
      }
      dci_.allocate(nMonomer-1);
      eta_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         dci_[i].allocate(meshDimensions);
         eta_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);
      dwp_.allocate(meshDimensions);
   }

   template <int D>
   void PredCorrBdStep<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(wp_.capacity() == nMonomer);
      UTIL_CHECK(wf_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(wp_[i].capacity() == meshSize);
         UTIL_CHECK(wf_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dci_.capacity() == nMonomer-1);
      UTIL_CHECK(eta_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dci_[i].capacity() == meshSize);
         UTIL_CHECK(eta_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      UTIL_CHECK(dwp_.capacity() == meshSize);
   }

   template <int D>
   bool PredCorrBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Save current state
      simulator().saveState();

      // Copy current W fields from parent system
      for (i = 0; i < nMonomer; ++i) {
         wp_[i] = system().w().rgrid(i);
         wf_[i] = wp_[i];
      }

      // Store initial value of pressure field at all grid points in dwp_
      dwp_ = simulator().wc(nMonomer-1);

      // Define constants used for step
      const double vSystem = system().domain().unitCell().volume();
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);

      // Construct all random displacement (noise) components
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> & eta = eta_[j];
         for (k = 0; k < meshSize; ++k) {
            eta[k] = b*random().gaussian();
         }
      }

      // Compute predicted state wp_, and store initial force dci_

      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         RField<D> const & eta = eta_[j];
         RField<D> & dci = dci_[j];
         for (k = 0; k < meshSize; ++k) {
            dwc_[k] = a*dc[k] + eta[k];
            dci[k] = dc[k];
         }
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wp = wp_[i];
            evec = simulator().chiEvecs(j,i);
            for (k = 0; k < meshSize; ++k) {
               wp[k] += evec*dwc_[k];
            }
         }
      }

      // Set modified system fields at predicted state wp_
      system().setWRGrid(wp_);

      // Set function return value to indicate failure by default
      bool isConverged = false;

      // Enforce incompressibility at predicted state 
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else {
         UTIL_CHECK(system().hasCFields());

         // Compute components and derivatives at wp_
         simulator().clearData();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();

         // Compute change dwp_ in pressure field 
         // Note: On entry, dwp_ is the old pressure field
         RField<D> const & wp = simulator().wc(nMonomer-1);
         for (k = 0; k < meshSize; ++k) {
            dwp_[k] = wp[k] - dwp_[k];
         }

         // Adjust predicted pressure field to final monomer fields
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wf = wf_[i];
            for (k = 0; k < meshSize; ++k) {
               wf[k] += dwp_[k];
            }
         }

         // Full step (corrector) change in exchange fields
         const double ha = 0.5*a;
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D> const & dcp = simulator().dc(j);
            RField<D> const & dci = dci_[j];
            RField<D> const & eta = eta_[j];
            for (k = 0; k < meshSize; ++k) {
               dwc_[k] = ha*( dci[k] + dcp[k]) + eta[k];
            }
            for (i = 0; i < nMonomer; ++i) {
               RField<D> & wf = wf_[i];
               evec = simulator().chiEvecs(j,i);
               for (k = 0; k < meshSize; ++k) {
                  wf[k] += evec*dwc_[k];
               }
            }
         }

         // Set system fields after predictor step
         system().setWRGrid(wf_);

         // Apply compressor to final state
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

      // True iff compression was successful after predictor and corrector
      return isConverged;
   }

}
}
#endif
