#ifndef RP_CHI_DERIVATIVE_TPP
#define RP_CHI_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"
#include <pscf/interaction/Interaction.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   ChiDerivative<D,T>::ChiDerivative(typename T::Simulator& simulator,
                                     typename T::System& system)
    : AverageAnalyzerT(simulator, system)
   {  ParamComposite::setClassName("ChiDerivative"); }


   /*
   * Compute and return derivative w/ respect to chi.
   */
   template <int D, class T>
   double ChiDerivative<D,T>::compute()
   {
      // Preconditions
      UTIL_CHECK(system().w().hasData());
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);
      // Only valid for nMonomer == 2

      const double vSystem = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();
      double chi = system().interaction().chi(0,1);

      // Pre-requisite computations
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      // Get field Hamiltonian per monomer
      double hField = simulator().fieldHamiltonian();
      hField /= (double)nMonomerSystem;

      // Compute derivative of f (per monomer) w/respect to bare chi
      double dfdchi = -(hField - 0.5*simulator().sc(nMonomer - 1))/chi
                    + 0.25;

      // Convert to derivative of total free energy F
      dfdchi *= (double)nMonomerSystem;

      // With N term
      dfdchi += double(meshSize)/(2.0 * chi);

      return dfdchi;
   }

   #if 0
   template <int D, class T>
   void ChiDerivative<D,T>::outputValue(int step, double value)
   {
      int nSamplePerOutput = AverageAnalyzerT::nSamplePerOutput();
      if (simulator().hasRamp() && nSamplePerOutput == 1) {
         std::ofstream& outputFile = AverageAnalyzerT::outputFile_;
         UTIL_CHECK(outputFile.is_open());
         double chi = system().interaction().chi(0,1);
         outputFile << Int(step);
         outputFile << Dbl(chi);
         outputFile << Dbl(value) << "\n";
      } else {
         AverageAnalyzerT::outputValue(step, value);
      }
   }
   #endif

}
}
#endif
