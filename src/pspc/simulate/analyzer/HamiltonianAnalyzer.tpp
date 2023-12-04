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
#include <pspc/simulate/McSimulator.h>
//#include <util/accumulators/Average.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(McSimulator<D>& mcSimulator, 
                                               System<D>& system)
    : AverageListAnalyzer<D>(system),
      mcSimulatorPtr_(&mcSimulator),
      systemPtr_(&(mcSimulator.system())),
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
      if (!mcSimulator().hasWc()){
         if (!hasAnalyzeChi_){
            mcSimulator().analyzeChi();
            hasAnalyzeChi_ = true;
         }
         system().compute();
         mcSimulator().computeWc();
         mcSimulator().computeHamiltonian();
      }
      double ideal = mcSimulator().idealHamiltonian();
      // outputFile_ << Dbl(ideal, 20)
      setValue(idealId_, ideal);
   
      double field = mcSimulator().fieldHamiltonian();
      // outputFile_ << Dbl(field, 20)
      setValue(fieldId_, field);
   
      double total = mcSimulator().hamiltonian();
      // outputFile_ << Dbl(total, 20) << std::endl;
      setValue(totalId_, total);
   }
   
}
}
#endif
