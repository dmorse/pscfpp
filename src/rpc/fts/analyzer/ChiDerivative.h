#ifndef RPC_CHI_DERIVATIVE_H
#define RPC_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpc
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * \see \ref rpc_ChiDerivative_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
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

   };

   // Explicit instantiation declarations
   extern template class ChiDerivative<1>;
   extern template class ChiDerivative<2>;
   extern template class ChiDerivative<3>;

}
}
#endif
