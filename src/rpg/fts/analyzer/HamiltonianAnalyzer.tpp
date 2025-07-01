#ifndef RPG_HAMILTONIAN_ANALYZER_TPP
#define RPG_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"
#include <rpg/System.h>
#include <rpg/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpg
{
   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(Simulator<D>& simulator, System<D>& system)
    : AverageListAnalyzer<D>(system),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      hasAnalyzeChi_(false),
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

      idealId_ = 0;
      fieldId_ = 1;
      totalId_ = 2;
      AverageListAnalyzer<D>::initializeAccumulators(3);
 
      setName(idealId_, "ideal");
      setName(fieldId_, "field");
      setName(totalId_, "total");
   }

   /*
   * Output energy to file
   */
   template <int D>
   void HamiltonianAnalyzer<D>::compute() 
   {  
      UTIL_CHECK(system().w().hasData());
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      
      #if 0
      if (!simulator().hasWc()){
         if (!hasAnalyzeChi_){
            simulator().analyzeChi();
            hasAnalyzeChi_ = true;
         }
         system().compute();
         simulator().computeWc();
      }
      #endif

      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }

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
