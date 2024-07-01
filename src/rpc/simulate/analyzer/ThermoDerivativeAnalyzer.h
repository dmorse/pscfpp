#ifndef RPC_THERMO_DERIVATIVE_ANALYZER_H
#define RPC_THERMO_DERIVATIVE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/accumulators/Average.h>           // member

namespace Pscf {
namespace Rpc 
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Base class for analyzers that evaluate derivatives of Hamiltonian.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class ThermoDerivativeAnalyzer : public Analyzer<D>
   {
   
   public:
   
      /**
      * Constructor.
      */
      ThermoDerivativeAnalyzer(Simulator<D>& simulator, 
                               System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~ThermoDerivativeAnalyzer(); 

      /**
      * Read parameters from archive.
      * 
      * Default implementation, reads interval, (optionally) hasAverage, 
      * and (optionally) outputFileName.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
      
      /**
      * Setup before simulation loop.
      */
      virtual void setup();
      
      /**
      * Compute and return the thermodynamic derivative.
      */
      virtual double computeDerivative() = 0;
      
      /**
      * Return current variable value
      */
      virtual double variable() = 0;
      
      /**
      * Return derivative parameter type.
      */
      virtual std::string parameterType() = 0;
   
      /**
      * Compute a sampled value and update the accumulator.
      */
      virtual void sample(long iStep);

      /**
      * Write final results to file after a simulation.
      */
      virtual void output();
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      
   protected:
        
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /**
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();
      
      // Output file stream
      std::ofstream outputFile_;
      
      // Output filename
      std::string outputFileName_;
      
      // Statistical accumulator.
      Average accumulator_;
      
   private:
            
      /// Whether the Average object is needed?
      bool hasAverage_;
      
      /**
      * Does the parameter file have the outputFileName parameter?
      * 
      * If true, the analyzer would output every derivative value to the file.
      */ 
      bool hasOutputFile_;
      
      /// Pointer to parent Simulator
      Simulator<D>* simulatorPtr_;
      
      /// Pointer to the parent system.
      System<D>* systemPtr_;  
   };
   
   // Inline functions
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& ThermoDerivativeAnalyzer<D>::simulator()
   {  return *simulatorPtr_; }
   
   // Get the parent system.
   template <int D>
   inline System<D>& ThermoDerivativeAnalyzer<D>::system()
   {  return *systemPtr_; }
   
   // Suppress implicit instantiation
   #ifndef RPC_THERMO_DERIVATIVE_ANALYZER_TPP
   extern template class ThermoDerivativeAnalyzer<1>;
   extern template class ThermoDerivativeAnalyzer<2>;
   extern template class ThermoDerivativeAnalyzer<3>;
   #endif

}
}
#endif 
