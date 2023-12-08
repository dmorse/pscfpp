#ifndef PSPC_HAMILTONIAN_ANALYZER_TPP
#define PSPC_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"

#include <pspc/System.h>
#include <pspc/simulate/Simulator.h>
//#include <util/accumulators/Average.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(Simulator<D>& simulator, 
                                               System<D>& system)
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
      if (!simulator().hasWc()){
         if (!hasAnalyzeChi_){
            simulator().analyzeChi();
            hasAnalyzeChi_ = true;
         }
         system().compute();
         simulator().computeWc();
         simulator().computeHamiltonian();
      }
      double ideal = simulator().idealHamiltonian();
      // outputFile_ << Dbl(ideal, 20)
      setValue(idealId_, ideal);
   
      double field = simulator().fieldHamiltonian();
      // outputFile_ << Dbl(field, 20)
      setValue(fieldId_, field);
   
      double total = simulator().hamiltonian();
      // outputFile_ << Dbl(total, 20) << std::endl;
      setValue(totalId_, total);
   }
   
}
}
#endif
