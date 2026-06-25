#ifndef RP_EXPLICIT_BD_STEP_TPP
#define RP_EXPLICIT_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExplicitBdStep.h"

#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/system/System.h>
//#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
//#include <rp/field/WFields.h>

#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   ExplicitBdStep<D,T>::ExplicitBdStep(BdSimulator<D,T>& simulator)
    : BdStepT(simulator),
      w_(),
      dwc_(),
      gaussianField_(),
      mobility_(0.0)
   {  ParamComposite::setClassName("ExplicitBdStep"); }

   /*
   * Read body of parameter file block and allocate memory.
   */
   template <int D, class T>
   void ExplicitBdStep<D,T>::readParameters(std::istream &in)
   {
      ParamComposite::read(in, "mobility", mobility_);

      // Allocate memory
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);
   }

   /*
   * Setup before entering simulation loop.
   */
   template <int D, class T>
   void ExplicitBdStep<D,T>::setup()
   {
      // Check array capacities
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      gaussianField_.allocate(meshDimensions);
   }

   /*
   * Take a single BD step.
   */
   template <int D, class T>
   bool ExplicitBdStep<D,T>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;

      // Save current state
      simulator().saveState();

      // Copy current W fields from parent system into w_
      for (i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);

      // Constants for normal distribution
      const double stddev = 1.0;
      const double mean = 0.0;

      // Modify local field copy wc_
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {

         // Generate normal distributed random numbers
         vecRandom().normal(gaussianField_, stddev, mean);

         // Compute change dwc_
         typename T::RField const & dc = simulator().dc(j);
         VecOp::addVcVc(dwc_, dc, a, gaussianField_, b);

         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            typename T::RField & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(w, dwc_, evec);
         }

      }

      // Set modified fields in parent system
      system().w().setRGrid(w_);
      simulator().clearData();

      // Enforce incompressibility (also solves MDE repeatedly)
      bool isConverged = false;
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else {
         isConverged = true;
         UTIL_CHECK(system().c().hasData());

         // Evaluate component properties in new state
         simulator().clearState();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();
      }

      return isConverged;
   }

}
}
#endif
