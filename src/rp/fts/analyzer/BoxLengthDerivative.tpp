#ifndef RP_BOX_LENGTH_DERIVATIVE_TPP
#define RP_BOX_LENGTH_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BoxLengthDerivative.h"
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   BoxLengthDerivative<D,T>::BoxLengthDerivative(
                                      typename T::Simulator& simulator,
                                      typename T::System& system)
    : AverageAnalyzerT(simulator, system)
   {  ParamComposite::setClassName("BoxLengthDerivative"); }

   /*
   * Compute and return required derivative.
   */
   template <int D, class T>
   double BoxLengthDerivative<D,T>::compute()
   {
      UTIL_CHECK(system().w().hasData());

      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      // Simulations for cubic box
      UTIL_CHECK(D == 3);

      int nParameter = system().domain().unitCell().nParameter();
      UTIL_CHECK(nParameter == 1);

      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();

      // Compute field Hamiltonian per monomer.
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      // Box length
      double l = system().domain().unitCell().parameter(0);

      // Obtain fieldHamiltonian
      double HW = simulator().fieldHamiltonian();

      // The fieldHamiltonian contribution to the derivative
      double dFdL = HW/vSystem * 3.0 * l * l;

      // Obtain ideal gas Hamiltonian - n[lnQ_id + W_{+}/M]
      double nlnQ= simulator().idealHamiltonian();
      dFdL += 3.0/l * nlnQ;

      // Obtain stress -1/Q dQdl per monomer
      if (!system().mixture().hasStress()) {
         system().computeStress();
      }
      double stress = system().mixture().stress(0);
      dFdL += stress * nMonomerSystem;

      // With N term
      dFdL -= 3.0 * double(meshSize)/(2.0 * l);

      return dFdL;
   }

   template <int D, class T>
   void BoxLengthDerivative<D,T>::outputValue(int step, double value)
   {
      if (simulator().hasRamp()
          && AverageAnalyzerT::nSamplePerOutput() == 1)
      {
         double l = system().domain().unitCell().parameter(0);

         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(l);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
      } else {
         AverageAnalyzerT::outputValue(step, value);
      }
   }

}
}
#endif
