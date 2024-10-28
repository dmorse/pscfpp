#ifndef RPC_PERTURBATION_DERIVATIVE_H
#define RPC_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ThermoDerivativeAnalyzer.h"
#include <rpc/System.h>
#include <rpc/fts/Simulator.h>

namespace Pscf {
namespace Rpc 
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H w/ respect to perturbation parameter lambda.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class PerturbationDerivative : public ThermoDerivativeAnalyzer<D>
   {
   
   public:
   
      /**
      * Constructor.
      */
      PerturbationDerivative(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~PerturbationDerivative();
   
      /**
      * Read parameters from archive.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
      
      /**
      * Setup before simulation loop.
      */
      virtual void setup();
      
      /**
      * Compute and return the derivative of H w/ respect to lambda.
      */
      virtual double computeDerivative();
      
      /**
      * Return current lambda value.
      */
      virtual double variable();
      
      /**
      * Return derivative parameter type "Perturbation Derivatie".
      */
      virtual std::string parameterType();
      
      using ParamComposite::setClassName;
      
   protected:
 
      using ThermoDerivativeAnalyzer<D>::simulator;
      using ThermoDerivativeAnalyzer<D>::system;         
   };
   
   // Suppress implicit instantiation
   #ifndef RPC_PERTURBATION_DERIVATIVE_TPP
   extern template class PerturbationDerivative<1>;
   extern template class PerturbationDerivative<2>;
   extern template class PerturbationDerivative<3>;
   #endif

}
}
#endif 
