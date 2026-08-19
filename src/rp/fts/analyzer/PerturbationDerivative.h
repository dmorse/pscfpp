#ifndef RP_PERTURBATION_DERIVATIVE_H
#define RP_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h>  // base class template
#include <pscf/backends/TmplDeclare.h>

namespace Pscf {
namespace Rp {

   /**
   * Evaluate derivative of H w/ respect to perturbation parameter lambda.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see rp_PerturbationDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class PerturbationDerivative : public AverageAnalyzer<D,T>
   {

   public:

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

   protected:

      /**
      * Compute and return the derivative of H w/ respect to lambda.
      */
      double compute() override;

      // Aliases for base classes.
      using AnalyzerT = Analyzer<D,T>;
      using AverageAnalyzerT = AverageAnalyzer<D,T>;

      // Inherited protected member functions (selected).
      using Analyzer<D,T>::simulator;
      using Analyzer<D,T>::system;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(PerturbationDerivative)

} // namespace Rp
} // namespace Pscf
#endif
