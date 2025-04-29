#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
      
      using AverageAnalyzer<D>::readParameters;   
      using AverageAnalyzer<D>::nSamplePerOutput;    
      using AverageAnalyzer<D>::setup;  
      using AverageAnalyzer<D>::sample;         
      using AverageAnalyzer<D>::output;         
      
   protected:
 
      using ParamComposite::setClassName;     
      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;
      using AverageAnalyzer<D>::outputFile_;
    
      /**
      * Compute and return the derivative of H w/ respect to concentration.
      */
      virtual double compute();
      
      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);
      
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
