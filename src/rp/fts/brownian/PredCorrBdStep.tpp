#ifndef RP_PRED_CORR_BD_STEP_TPP
#define RP_PRED_CORR_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   PredCorrBdStep<D,T>::PredCorrBdStep(BdSimulator<D,T>& simulator)
    : BdStepT(simulator),
      wf_(),
      dci_(),
      eta_(),
      dwc_(),
      dwp_(),
      mobility_(0.0)
   {  ParamComposite::setClassName("PredCorrBdStep"); }

   /*
   * Read body of parameter file block and allocate memory.
   */
   template <int D, class T>
   void PredCorrBdStep<D,T>::readParameters(std::istream &in)
   {
      ParamComposite::read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      wp_.allocate(nMonomer);
      wf_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         wp_[i].allocate(meshDimensions);
         wf_[i].allocate(meshDimensions);
      }
      dci_.allocate(nMonomer-1);
      eta_.allocate(nMonomer-1);
      for (int i = 0; i < nMonomer - 1; ++i) {
         dci_[i].allocate(meshDimensions);
         eta_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);
      dwp_.allocate(meshDimensions);
   }

   /*
   * Setup before entering main simulation loop.
   */
   template <int D, class T>
   void PredCorrBdStep<D,T>::setup()
   {
      // Check array capacities
      const int meshSize = system().domain().mesh().size();
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(wp_.capacity() == nMonomer);
      UTIL_CHECK(wf_.capacity() == nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(wp_[i].capacity() == meshSize);
         UTIL_CHECK(wf_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dci_.capacity() == nMonomer - 1);
      UTIL_CHECK(eta_.capacity() == nMonomer - 1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dci_[i].capacity() == meshSize);
         UTIL_CHECK(eta_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      UTIL_CHECK(dwp_.capacity() == meshSize);
   }

   /*
   * Take one BD step.
   */
   template <int D, class T>
   bool PredCorrBdStep<D,T>::step()
   {
      // Save initial state
      simulator().saveState();

      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;

      // Copy initial w monomer fields into wp_ and wf_
      for (i = 0; i < nMonomer; ++i) {
         wp_[i] = system().w().rgrid(i);
         wf_[i] = wp_[i];
      }

      // Copy initial value of pressure field as dwp_
      dwp_ = simulator().wc(nMonomer-1);

      // Constants used for step
      const double vSystem = system().domain().unitCell().volume();
      const double a = -1.0 * mobility_;
      const double stddev = sqrt(2.0*mobility_*double(meshSize)/vSystem);
      const double mean = 0.0;

      // Generate all random displacement components
      for (j = 0; j < nMonomer - 1; ++j) {
         BdStepT::vecRandom().normal(eta_[j], stddev, mean);
      }

      // Compute predicted state wp_, and store initial force as dci_

      // Loop over composition eigenvectors of projected chi matrix
      for (j = 0; j < nMonomer - 1; ++j) {
         RFieldT const & dc = simulator().dc(j);
         VecOp::addVcVc(dwc_, dc, a, eta_[j], 1.0);
         VecOp::eqV(dci_[j], dc);
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            double evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(wp_[i], dwc_, evec);
         }
      }

      // Set system fields to predicted state wp_
      system().w().setRGrid(wp_);

      // Apply compressor after predictor step
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
         bool isConverged = false;
         return isConverged;
      }

      // Compute components and derivatives at wp_
      UTIL_CHECK(system().c().hasData());
      simulator().clearData();
      simulator().computeWc();
      simulator().computeCc();
      simulator().computeDc();

      // Compute change dwp_ in pressure field
      // Note: On entry, dwp_ is the old pressure field
      RFieldT const & wp = simulator().wc(nMonomer-1);
      VecOp::subVV(dwp_, wp, dwp_);

      // Adjust predicted pressure field
      for (i = 0; i < nMonomer; ++i) {
         VecOp::addEqV(wf_[i], dwp_);
      }

      // Compute corrected state wf_
      const double ha = 0.5*a;
      for (j = 0; j < nMonomer - 1; ++j) {
         RFieldT const & dcp = simulator().dc(j);
         VecOp::addVcVcVc(dwc_, dci_[j], ha, dcp, ha, eta_[j], 1.0);
         for (i = 0; i < nMonomer; ++i) {
            double evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(wf_[i], dwc_, evec);
         }
      }

      // Set system fields after corrector step
      system().w().setRGrid(wf_);

      // Apply compressor to final state
      int compress2 = simulator().compressor().compress();
      if (compress2 != 0){
         simulator().restoreState();
         bool isConverged = false;
         return isConverged;
      }

      // Compute components and derivatives in final state
      UTIL_CHECK(system().c().hasData());
      simulator().clearState();
      simulator().clearData();
      simulator().computeWc();
      simulator().computeCc();
      simulator().computeDc();

      // Success
      bool isConverged = true;
      return isConverged;
   }

}
}
#endif
