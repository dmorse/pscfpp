#ifndef RP_HAMILTONIAN_ANALYZER_TPP
#define RP_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   HamiltonianAnalyzer<D,T>::HamiltonianAnalyzer(
                                      Simulator<D,T>& simulator,
                                      System<D,T>& system)
    : AverageListAnalyzer<D,T>(simulator, system),
      idealId_(-1),
      fieldId_(-1),
      totalId_(-1)
   {  ParamComposite::setClassName("HamiltonianAnalyzer"); }

   /*
   * Read interval and outputFileName.
   */
   template <int D, class T>
   void HamiltonianAnalyzer<D,T>::readParameters(std::istream& in)
   {
      AverageListAnalyzer<D,T>::readParameters(in);
      AverageListAnalyzer<D,T>::initializeAccumulators(3);

      idealId_ = 0;
      AverageListAnalyzer<D,T>::setName(idealId_, "ideal");
      fieldId_ = 1;
      AverageListAnalyzer<D,T>::setName(fieldId_, "field");
      totalId_ = 2;
      AverageListAnalyzer<D,T>::setName(totalId_, "total");
   }

   /*
   * Output energy to file.
   */
   template <int D, class T>
   void HamiltonianAnalyzer<D,T>::compute()
   {
      UTIL_CHECK(system().w().hasData());
      if (!system().c().hasData()) {
         system().compute();
      }
      UTIL_CHECK(system().c().hasData());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      UTIL_CHECK(simulator().hasWc());
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      UTIL_CHECK(simulator().hasHamiltonian());

      double ideal = simulator().idealHamiltonian();
      AverageListAnalyzer<D,T>::setValue(idealId_, ideal);

      double field = simulator().fieldHamiltonian();
      AverageListAnalyzer<D,T>::setValue(fieldId_, field);

      double total = simulator().hamiltonian();
      AverageListAnalyzer<D,T>::setValue(totalId_, total);
   }

}
}
#endif
