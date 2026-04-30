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

} // namespace Rp
} // namespace Pscf
#endif
