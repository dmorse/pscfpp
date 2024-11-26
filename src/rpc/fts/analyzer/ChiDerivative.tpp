#ifndef RPC_CHI_DERIVATIVE_TPP
#define RPC_CHI_DERIVATIVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpc 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ChiDerivative<D>::ChiDerivative(Simulator<D>& simulator, 
                                                     System<D>& system) 
    : ThermoDerivativeAnalyzer<D>(simulator, system)
   { setClassName("ChiDerivative"); }

   /*
   * Destructor.
   */
   template <int D>
   ChiDerivative<D>::~ChiDerivative() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void ChiDerivative<D>::readParameters(std::istream& in) 
   {
      ThermoDerivativeAnalyzer<D>::readParameters(in);
   }
   
   /*
   * Setup before simulation loop.
   */ 
   template <int D>
   void ChiDerivative<D>::setup()
   {}
   
   template <int D>
   double ChiDerivative<D>::computeDerivative()
   {
      UTIL_CHECK(system().w().hasData());
      
      // For AB diblock
      const int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer == 2);
      
      const double vSystem  = system().domain().unitCell().volume();
      const double vMonomer = system().mixture().vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      const int meshSize = system().domain().mesh().size();
      double chi= system().interaction().chi(0,1);
      
      /* 
      * Obteain field Hamiltonian per monomer.
      * The fieldHamitonian is calculated in the computeHamiltonian() function,
      * located in rpc/fts/Simulator.tpp 
      */
      if (!system().hasCFields()) {
         system().compute();
      }
      if (!simulator().hasWc()){
         simulator().computeWc();
      }
      if (!simulator().hasHamiltonian()) {
         simulator().computeHamiltonian();
      }
      double hField = simulator().fieldHamiltonian()/nMonomerSystem;
      
      // Compute derivative of free energy with repect to chi_bare per monomer
      double dfdchi = -(hField - 0.5*simulator().sc(nMonomer - 1))/chi + 1.0/4.0;
      
      dfdchi *= nMonomerSystem;
 
      // With N term
      dfdchi += double(meshSize)/(2.0 * chi);
      
      return dfdchi;
   }
   
   template <int D>
   double ChiDerivative<D>::variable()
   { return system().interaction().chi(0,1);  }
   
   template <int D>
   std::string ChiDerivative<D>::parameterType()
   { return "Chi Derivative";  }

}
}
#endif
