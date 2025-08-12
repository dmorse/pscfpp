#ifndef RPG_CHI_DERIVATIVE_H
#define RPG_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpg 
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class ChiDerivative : public AverageAnalyzer<D>
   {
   
   public:
   
      /**
      * Constructor.
      */
      ChiDerivative(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~ChiDerivative(); 

      /**
      * Compute and return the derivative of H w/ respect to chi.
      */
      virtual double compute();
      
      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);
      
      using AverageAnalyzer<D>::readParameters;
      using AverageAnalyzer<D>::nSamplePerOutput;
      using AverageAnalyzer<D>::setup;
      using AverageAnalyzer<D>::sample;
      using AverageAnalyzer<D>::output; 
      
   protected:
 
      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system; 
      using AverageAnalyzer<D>::outputFile_;
      using ParamComposite::setClassName;

   };
   
   #ifndef RPG_CHI_DERIVATIVE_TPP
   // Suppress implicit instantiation
   extern template class ChiDerivative<1>;
   extern template class ChiDerivative<2>;
   extern template class ChiDerivative<3>;
   #endif

}
}
#endif 
