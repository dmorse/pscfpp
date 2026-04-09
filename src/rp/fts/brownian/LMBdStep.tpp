#ifndef RP_LM_BD_STEP_TPP
#define RP_LM_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LMBdStep.h"

#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   LMBdStep<D,T>::LMBdStep(typename T::BdSimulator& simulator)
    : BdStepT(simulator),
      w_(),
      etaA_(),
      etaB_(),
      dwc_(),
      etaNewPtr_(nullptr),
      etaOldPtr_(nullptr),
      mobility_(0.0)
   {  ParamComposite::setClassName("LmBdStep"); }

   /*
   * Read body of parameter file block and allocate memory.
   */
   template <int D, class T>
   void LMBdStep<D,T>::readParameters(std::istream &in)
   {
      ParamComposite::read(in, "mobility", mobility_);

      // Allocate memory
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
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
   }

   /*
   * Generate new random displacement values.
   */
   template <int D, class T>
   void LMBdStep<D,T>::generateEtaNew()
   {
      // Constants
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vSystem = system().domain().unitCell().volume();
      const double stddev = sqrt(0.5*mobility_*double(meshSize)/vSystem);
      const double mean = 0.0;

      for (int j = 0; j < nMonomer - 1; ++j) {
         vecRandom().normal(etaNew(j), stddev, mean);
      }
   }

   /*
   * Exchange pointers to old and new random fields.
   */
   template <int D, class T>
   void LMBdStep<D,T>::exchangeOldNew()
   {
      DArray< RFieldT >* temp;
      temp = etaOldPtr_;
      etaOldPtr_ = etaNewPtr_;
      etaNewPtr_ = temp;
   }

   /*
   * Setup before main simulation loop.
   */
   template <int D, class T>
   void LMBdStep<D,T>::setup()
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
   template <int D, class T>
   bool LMBdStep<D,T>::step()
   {
      // Save current state
      simulator().saveState();

      // Copy current W fields from parent system into w_
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Generate new random displacement values
      generateEtaNew();

      // Take LM step
      const double a = -1.0*mobility_;
      double evec;
      // Loop over eigenvectors of projected chi matrix
      for (int j = 0; j < nMonomer - 1; ++j) {
         RFieldT const & etaN = etaNew(j);
         RFieldT const & etaO = etaOld(j);
         RFieldT const & dc = simulator().dc(j);
         VecOp::addVV(dwc_, etaN, etaO);
         VecOp::addEqVc(dwc_, dc, a);
         // Loop over monomer types
         for (int i = 0; i < nMonomer; ++i) {
            evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(w_[i], dwc_, evec);
         }
      }

      // Set modified fields in parent system
      system().w().setRGrid(w_);

      // Enforce incompressibility (also solves MDE repeatedly)
      bool isConverged = false;
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else {
         isConverged = true;
         UTIL_CHECK(system().c().hasData());

         // Compute components and derivatives at wp_
         simulator().clearState();
         simulator().clearData();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();

         // Exchange old and new random fields
         exchangeOldNew();
      }

      return isConverged;
   }

}
}
#endif
