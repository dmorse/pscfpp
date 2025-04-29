#ifndef RPC_HAMILTONIAN_ANALYZER_TPP
#define RPC_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>

#include <iostream>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(Simulator<D>& simulator, 
                                               System<D>& system)
    : AverageListAnalyzer<D>(simulator, system),
      idealId_(-1),
      fieldId_(-1),
      totalId_(-1)
   {  setClassName("HamiltonianAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void HamiltonianAnalyzer<D>::readParameters(std::istream& in) 
   {
      AverageListAnalyzer<D>::readParameters(in);

      AverageListAnalyzer<D>::initializeAccumulators(3);
 
      idealId_ = 0;
      setName(idealId_, "ideal");
      fieldId_ = 1;
      setName(fieldId_, "field");
      totalId_ = 2;
      setName(totalId_, "total");
   }

   /*
   * Output energy to file
   */
   template <int D>
   void HamiltonianAnalyzer<D>::compute() 
   {
      UTIL_CHECK(system().w().hasData());

      if (!system().hasCFields()) {
         system().compute();
      }
      UTIL_CHECK(system().hasCFields());
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      UTIL_CHECK(simulator().hasWc());
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      UTIL_CHECK(simulator().hasHamiltonian());

      double ideal = simulator().idealHamiltonian();
      setValue(idealId_, ideal);
   
      double field = simulator().fieldHamiltonian();
      setValue(fieldId_, field);
   
      double total = simulator().hamiltonian();
      setValue(totalId_, total);
   }
   
}
}
#endif
