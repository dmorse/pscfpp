#ifndef RP_PERTURBATION_DERIVATIVE_H
#define RP_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h>  // base class template

namespace Pscf {
namespace Rp {

   /**
   * Evaluate derivative of H w/ respect to perturbation parameter lambda.
   *
   * Specializations of this template are used as base classes for two
   * closely analogous class templates, also named PerturbationDerivative, 
   * that are defined in the Rpc and Rpg namespaces for use in the 
   * pscf_rpc and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see rp_PerturbationDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class PerturbationDerivative : public T::AverageAnalyzer
   {

   protected:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      PerturbationDerivative(Simulator<D,T>& simulator, 
                             System<D,T>& system);

      /**
      * Destructor.
      */
      ~PerturbationDerivative() = default;

      /**
      * Compute and return the derivative of H w/ respect to lambda.
      */
      double compute() override;

      // Aliases for base classes.
      using AnalyzerT = typename T::Analyzer;
      using AverageAnalyzerT = typename T::AverageAnalyzer;

      // Inherited protected member functions (selected).
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   };

} // namespace Rp
} // namespace Pscf
#endif
