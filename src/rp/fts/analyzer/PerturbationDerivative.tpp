#ifndef RP_PERTURBATION_DERIVATIVE_TPP
#define RP_PERTURBATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   PerturbationDerivative<D,T>::PerturbationDerivative(
                                        typename T::Simulator& simulator,
                                        typename T::System& system)
    : AverageAnalyzerT(simulator, system)
   {  ParamComposite::setClassName("PerturbationDerivative"); }

   /*
   * Compute and return derivative df.
   */
   template <int D, class T>
   double PerturbationDerivative<D,T>::compute()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(simulator().hasPerturbation());

      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      return simulator().perturbation().df();
   }

   #if 0
   template <int D, class T>
   void PerturbationDerivative<D,T>::outputValue(int step, double value)
   {
      int nSamplePerOutput = AverageAnalyzerT::nSamplePerOutput();
      if (simulator().hasRamp() && nSamplePerOutput == 1) {
         std::ofstream& outputFile = AverageAnalyzerT::outputFile_;
         UTIL_CHECK(outputFile.is_open());
         double lambda = simulator().perturbation().lambda();
         outputFile << Int(step);
         outputFile << Dbl(lambda);
         outputFile << Dbl(value);
         outputFile << "\n";
      } else {
         AverageAnalyzerT::outputValue(step, value);
      }
   }
   #endif

} // namespace Rp
} // namespace Pscf
#endif
