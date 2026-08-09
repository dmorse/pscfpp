#ifndef RP_CONCENTRATION_DERIVATIVE_H
#define RP_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h> // base class template
#include <pscf/backends/TmplDeclare.h>

namespace Pscf {
namespace Rp {

   /**
   * Evaluate the derivative of H with respect to monomer concentration.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, CppTp<D> or CudaTp<D>.
   *
   * \see \ref rp_ConcentrationDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class ConcentrationDerivative : public AverageAnalyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationDerivative(Simulator<D,T>& simulator,
                              System<D,T>& system);

      /**
      * Destructor.
      */
      ~ConcentrationDerivative() = default;

   protected:

      /**
      * Compute and return the derivative of H w/ respect to concentration.
      */
      double compute() override;

      // Aliases for base classes.
      //using AnalyzerT = Analyzer<D,T>;
      //using AverageAnalyzerT = AverageAnalyzer<D,T>;

      // Inherited protected member functions (selected).
      using Analyzer<D,T>::simulator;
      using Analyzer<D,T>::system;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(ConcentrationDerivative)

} // namespace Rp
} // namespace Pscf
#endif
