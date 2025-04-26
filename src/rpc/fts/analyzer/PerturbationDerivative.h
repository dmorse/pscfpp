#ifndef RPC_PERTURBATION_DERIVATIVE_H
#define RPC_PERTURBATION_DERIVATIVE_H

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
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H w/ respect to perturbation parameter lambda.
   *
   * \see rpc_PerturbationDerivative_page "ParameterFileFormat"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class PerturbationDerivative : public AverageAnalyzer<D>
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
      * Compute and return the derivative of H w/ respect to lambda.
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

   // Suppress implicit instantiation
   #ifndef RPC_PERTURBATION_DERIVATIVE_TPP
   extern template class PerturbationDerivative<1>;
   extern template class PerturbationDerivative<2>;
   extern template class PerturbationDerivative<3>;
   #endif

}
}
#endif
