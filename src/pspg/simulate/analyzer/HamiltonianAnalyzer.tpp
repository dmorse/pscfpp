#ifndef PSPG_HAMILTONIAN_ANALYZER_TPP
#define PSPG_HAMILTONIAN_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAnalyzer.h"
#include <pspg/simulate/McSimulator.h>
#include <pspg/System.h>
#include <util/accumulators/Average.h>


namespace Pscf {
namespace Pspg
{
   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAnalyzer<D>::HamiltonianAnalyzer(McSimulator<D>& mcSimulator, System<D>& system)
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

      // Count number of values
      int id = 0;
      idealId_ = id;
      ++id; 
      fieldId_ = id;
      ++id;
      totalId_ = id;
      ++id; 
      AverageListAnalyzer<D>::initializeAccumulators(id);
 
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
      if (!mcSimulator().hasWC()){
         if (!hasAnalyzeChi_){
            mcSimulator().analyzeChi();
            hasAnalyzeChi_ = true;
         }
         system().compute();
         mcSimulator().computeWC();
         mcSimulator().computeHamiltonian();
      }
      double ideal = mcSimulator().idealHamiltonian();
      // outputFile_ << Dbl(ideal, 20)
      setValue(idealId_, ideal);
   
      double field = mcSimulator().mcFieldHamiltonian();
      // outputFile_ << Dbl(field, 20)
      setValue(fieldId_, field);
   
      double total = mcSimulator().hamiltonian();
      // outputFile_ << Dbl(total, 20) << std::endl;
      setValue(totalId_, total);
   }
   
}
}
#endif
