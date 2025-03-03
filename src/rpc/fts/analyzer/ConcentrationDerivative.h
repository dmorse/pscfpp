#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpc 
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to concentration.
   *
   * \see \ref rpc_ConcentrationDerivative_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationDerivative : public AverageAnalyzer<D>
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
      * Compute and return the derivative of H w/ respect to concentration.
      */
      virtual double compute();
      
      using ParamComposite::setClassName;
      using AverageAnalyzer<D>::setup;
      using AverageAnalyzer<D>::readParameters;         
      using AverageAnalyzer<D>::sample;         
      using AverageAnalyzer<D>::output;         
      
   protected:
 
      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;         
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
