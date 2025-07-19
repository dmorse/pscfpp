#ifndef RPC_CONCENTRATION_DERIVATIVE_TPP
#define RPC_CONCENTRATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/fts/perturbation/Perturbation.h>

namespace Pscf {
namespace Rpc 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ConcentrationDerivative<D>::ConcentrationDerivative(Simulator<D>& simulator, 
                                                     System<D>& system) 
    : AverageAnalyzer<D>(simulator, system)
   { setClassName("ConcentrationDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   ConcentrationDerivative<D>::~ConcentrationDerivative() 
   {}

   template <int D>
   double ConcentrationDerivative<D>::compute()
   { 
      UTIL_CHECK(system().w().hasData());
      
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2); 
      
      double vMonomer = system().mixture().vMonomer();
      const int meshSize = system().domain().mesh().size();
      
      if (!system().c().hasData()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      
      // Obtain Hamiltonian
      double h = simulator().hamiltonian();
      
      // Calculate derivative with respect to concentration
      double dfdc = h * vMonomer;
      
      // With N term
      double Hh = double(meshSize)/2.0 * vMonomer;
      dfdc -= Hh;
      
      return dfdc;
   }
   
   template <int D>
   void ConcentrationDerivative<D>::outputValue(int step, double value)
   {
      if (simulator().hasRamp() && nSamplePerOutput() == 1) {
         double vMonomer = system().mixture().vMonomer();
         
         UTIL_CHECK(outputFile_.is_open());
         outputFile_ << Int(step);
         outputFile_ << Dbl(vMonomer);
         outputFile_ << Dbl(value);
         outputFile_ << "\n";
       } else {
         AverageAnalyzer<D>::outputValue(step, value);
       }
   }
   
}
}
#endif
