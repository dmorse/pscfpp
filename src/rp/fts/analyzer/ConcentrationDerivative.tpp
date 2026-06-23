#ifndef RP_CONCENTRATION_DERIVATIVE_TPP
#define RP_CONCENTRATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   ConcentrationDerivative<D,T>::ConcentrationDerivative(
                                  Simulator<D,T>& simulator,
                                  System<D,T>& system)
    : AverageAnalyzerT(simulator, system)
   {  ParamComposite::setClassName("ConcentrationDerivative"); }

   /*
   * Compute and return the derivative of interest.
   */
   template <int D, class T>
   double ConcentrationDerivative<D,T>::compute()
   {
      UTIL_CHECK(system().w().hasData());

      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      double vMonomer = system().mixture().vMonomer();
      const int meshSize = system().domain().mesh().size();

      // Compute and retrieve Hamiltonian
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      double h = simulator().hamiltonian();

      // Calculate derivative with respect to concentration
      double dfdc = h * vMonomer;

      // With N term
      double Hh = double(meshSize)/2.0 * vMonomer;
      dfdc -= Hh;

      return dfdc;
   }

}
}
#endif
