#ifndef RPC_PERTURBATION_DERIVATIVE_TPP
#define RPC_PERTURBATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/perturbation/Perturbation.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   PerturbationDerivative<D>::PerturbationDerivative(Simulator<D>& simulator, 
                                                     System<D>& system) 
    : AverageAnalyzer<D>(simulator, system)
   {  setClassName("PerturbationDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   PerturbationDerivative<D>::~PerturbationDerivative() 
   {}

   template <int D>
   double PerturbationDerivative<D>::compute()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(simulator().hasPerturbation());

      if (!system().hasCFields()) {
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
   
   template <int D>
   void PerturbationDerivative<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double lambda = simulator().perturbation().lambda(); 
         
         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(lambda);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }
   
}
}
#endif
