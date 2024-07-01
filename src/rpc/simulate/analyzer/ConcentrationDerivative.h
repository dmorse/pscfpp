#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ThermoDerivativeAnalyzer.h"
#include <rpc/System.h>
#include <rpc/simulate/Simulator.h>

namespace Pscf {
namespace Rpc 
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to concentration.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class ConcentrationDerivative : public ThermoDerivativeAnalyzer<D>
   {
   
   public:
   
      /**
      * Constructor.
      */
      ConcentrationDerivative(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~ConcentrationDerivative(); 

      /**
      * Read parameters from archive.
      * 
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
      
      /**
      * Setup before simulation loop.
      */
      virtual void setup();
      
      /**
      * Compute and return the derivative of H w/ respect to concentration
      */
      virtual double computeDerivative();
      
      /**
      * Return current variable value
      */
      virtual double variable();
      
      /**
      * Return the derivative parameter type string "Concentration Derivative"
      */
      virtual std::string parameterType();
      
      using ParamComposite::setClassName;
      
   protected:
 
      using ThermoDerivativeAnalyzer<D>::simulator;
      using ThermoDerivativeAnalyzer<D>::system;         
   };
   
   // Suppress implicit instantiation
   #ifndef RPC_CONCENTRATION_DERIVATIVE_TPP
   extern template class ConcentrationDerivative<1>;
   extern template class ConcentrationDerivative<2>;
   extern template class ConcentrationDerivative<3>;
   #endif

}
}
#endif 
