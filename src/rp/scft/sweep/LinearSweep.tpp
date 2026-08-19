#ifndef RP_LINEAR_SWEEP_TPP
#define RP_LINEAR_SWEEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearSweep.h"
#include <rp/system/System.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   template <int D, class T>
   LinearSweep<D,T>::LinearSweep(System<D,T>& system)
    : SweepT(system)
   {}

   template <int D, class T>
   void LinearSweep<D,T>::readParameters(std::istream& in)
   {
      // Call the base class's readParameters function
      SweepT::readParameters(in);

      // Read the number of sweep parameters, allocate parameters_ array
      ParamComposite::read(in, "nParameter", nParameter_);
      parameters_.allocate(nParameter_);

      // Set the pointer to the array of specialized parameter types for
      // each SweepParameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setParameterTypesArray(SweepTmplT::parameterTypes_);
      }

      // Read array of SweepParameters, calling << operator for each
      ParamComposite::readDArray(in, "parameters", parameters_,nParameter_);

      // Verify net zero change in volume fractions, if these are swept
      double sum = 0.0;
      for (int i = 0; i < nParameter_; ++i) {
         if (parameters_[i].type() == "phi_polymer" ||
             parameters_[i].type() == "phi_solvent")
         {  sum += parameters_[i].change(); }
      }
      UTIL_CHECK(sum > -0.000001);
      UTIL_CHECK(sum < 0.000001);
   }

   template <int D, class T>
   void LinearSweep<D,T>::setup()
   {
      // Verify that the LinearSweep has a pointer to parent System
      UTIL_CHECK(SweepT::hasSystem());

      // Call base class setup function
      SweepT::setup();

      // Set system pointer and initial value for each parameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setSystem(system());
         parameters_[i].getInitial();
      }
   }

   template <int D, class T>
   void LinearSweep<D,T>::setParameters(double s)
   {
      // Update the system parameter values
      double newVal;
      for (int i = 0; i < nParameter_; ++i) {
         newVal = parameters_[i].initial() + s*parameters_[i].change();
         parameters_[i].update(newVal);
      }
   }

   template <int D, class T>
   void LinearSweep<D,T>::outputSummary(std::ostream& out)
   {}

}
}
#endif
