#ifndef RPG_CHI_DERIVATIVE_TPP
#define RPG_CHI_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"

#include <rpg/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <pscf/inter/Interaction.h>

namespace Pscf {
namespace Rpg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ChiDerivative<D>::ChiDerivative(Simulator<D>& simulator,
                                   System<D>& system)
    : AverageAnalyzer<D>(simulator, system)
   {  setClassName("ChiDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   ChiDerivative<D>::~ChiDerivative()
   {}

   template <int D>
   double ChiDerivative<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());

      // Check that nMonomer == 2
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);

      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();
      double chi = system().interaction().chi(0,1);

      // Pre-requisite computations
      if (!system().hasCFields()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

      // Get field Hamiltonian per monomer
      double hField = simulator().fieldHamiltonian()/nMonomerSystem;

      // Compute derivative of f (per monomer) w/respect to chi_bare
      double dfdchi = -(hField - 0.5*simulator().sc(nMonomer - 1))/chi
                    + 1.0/4.0;

      // Convert to derivative of total free energy F
      dfdchi *= nMonomerSystem;

      // With N term
      dfdchi += double(meshSize)/(2.0 * chi);

      return dfdchi;
   }
   
   template <int D>
   void ChiDerivative<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double chi= system().interaction().chi(0,1);
         
         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(chi);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
      } else {
         AverageAnalyzer<D>::outputValue(step, value);
      }
   }

}
}
#endif
