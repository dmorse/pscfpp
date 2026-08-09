#ifndef RP_CUBIC_LENGTH_DERIVATIVE_TPP
#define RP_CUBIC_LENGTH_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CubicLengthDerivative.h"

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   CubicLengthDerivative<D,T>::CubicLengthDerivative(
                                      Simulator<D,T>& simulator,
                                      System<D,T>& system)
    : AverageAnalyzer<D,T>(simulator, system)
   {  ParamComposite::setClassName("CubicLengthDerivative"); }

   /*
   * Compute and return required derivative.
   */
   template <int D, class T>
   double CubicLengthDerivative<D,T>::compute()
   {
      // Preconditions
      UTIL_CHECK(D == 3);
      UTIL_CHECK(system().domain().unitCell().nParameter() == 1);
      UTIL_CHECK(system().w().hasData());

      // Compute field Hamiltonian per monomer.
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!system().mixture().hasStress()) {
         system().computeStress();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      // Length of cubic box
      double l = system().domain().unitCell().parameter(0);

      // Term extensive in Hamiltonian
      double const hamiltonian = simulator().hamiltonian();
      double dFdL = hamiltonian * (3.0 / l);

      // "Stress" term arising from changes in ln Q
      double const stress = system().mixture().stress(0);
      double const vSystem  = system().domain().unitCell().volume();
      double const vMonomer = system().mixture().vMonomer();
      dFdL += stress * vSystem / vMonomer;

      // Term arising from ln N
      int const meshSize = system().domain().mesh().size();
      dFdL -= 3.0 * double(meshSize) / (2.0 * l);

      return dFdL;
   }

}
}
#endif
