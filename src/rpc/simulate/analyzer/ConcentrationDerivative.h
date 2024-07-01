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
   * Evaluate thermodynamic derivatives
   *
   * This class evaluates the derivative of free energy with respect
   * to relevant parameter.
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
      * Read interval (optionally) outputFileName 
      * and (optionally) hasAverage.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
      
      /**
      * Setup before loop. Opens an output file, if any.
      */
      virtual void setup();
      
      /**
      * Compute the thermodynamic derivative.
      */
      virtual double computeDerivative();
      
      /**
      * Obtain current variable value
      */
      virtual double variable();
      
      /**
      * Obtain derivative parameter type.
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
