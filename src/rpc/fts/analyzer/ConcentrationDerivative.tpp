#ifndef RPC_CONCENTRATION_DERIVATIVE_TPP
#define RPC_CONCENTRATION_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
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
    : ThermoDerivativeAnalyzer<D>(simulator, system)
   { setClassName("ConcentrationDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   ConcentrationDerivative<D>::~ConcentrationDerivative() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ConcentrationDerivative<D>::readParameters(std::istream& in) 
   {
      ThermoDerivativeAnalyzer<D>::readParameters(in);
   }
   
   /*
   * Setup before simulation loop.
   */ 
   template <int D>
   void ConcentrationDerivative<D>::setup()
   {}
   
   template <int D>
   double ConcentrationDerivative<D>::computeDerivative()
   { 
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2); 
      
      double vMonomer = system().mixture().vMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // Obteain field hamiltonian
      if (!simulator().hasWc()){
         system().compute();
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      
      // Obteain Hamiltonian
      double h = simulator().hamiltonian();
      
      // Calculate derivative with respect to concentration
      double dfdc = h * vMonomer;
      
      // With N term
      double Hh = double(meshSize)/2.0 * vMonomer;
      dfdc -= Hh;
      
      return dfdc;
   }
   
   template <int D>
   double ConcentrationDerivative<D>::variable()
   { return system().mixture().vMonomer(); }
   
   template <int D>
   std::string ConcentrationDerivative<D>::parameterType()
   { return "Concentration Derivative"; }

}
}
#endif
