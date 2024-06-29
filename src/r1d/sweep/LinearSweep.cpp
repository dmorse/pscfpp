/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearSweep.h"
#include <r1d/System.h>
#include <cstdio>

namespace Pscf {
namespace R1d {

   using namespace Util;

   LinearSweep::LinearSweep(System& system)
   : Sweep(system)
   {}

   void LinearSweep::readParameters(std::istream& in)
   {
      // Call the base class's readParameters function.
      Sweep::readParameters(in);
      
      // Read in the number of sweep parameters and allocate.
      read(in, "nParameter", nParameter_);
      parameters_.allocate(nParameter_);

      // Set the pointer to the array of specialized parameter types for 
      // each SweepParameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setParameterTypesArray(
                        SweepTmpl< DArray<System::WField> >::parameterTypes_);
      }
      
      // Read in array of SweepParameters, calling << for each
      readDArray< SweepParameter >(in, "parameters", 
                                   parameters_, nParameter_);

      // Verify net zero change in volume fractions if being swept
      double sum = 0.0;
      for (int i = 0; i < nParameter_; ++i) {
         if (parameters_[i].type() == "phi_polymer" 
             || parameters_[i].type() == "phi_solvent") 
         {
            sum += parameters_[i].change();
         }
      }
      UTIL_CHECK(sum > -0.000001);
      UTIL_CHECK(sum < 0.000001);
   }

   void LinearSweep::setup()
   {
      // Call base class's setup function
      Sweep::setup();
      
      // Set system pointer and initial value for each parameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setSystem(system());
         parameters_[i].getInitial();
      }
   }

   void LinearSweep::setParameters(double s)
   {
      // Update the system parameter values
      double newVal;
      for (int i = 0; i < nParameter_; ++i) {
         newVal = parameters_[i].initial() + s*parameters_[i].change();
         parameters_[i].update(newVal);
      }
   }

   void LinearSweep::outputSummary(std::ostream& out)
   {}

}
}

