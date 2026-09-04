#ifndef RP_CUBIC_LENGTH_DERIVATIVE_H
#define RP_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h> // base class template
#include <pscf/backend/TmplDeclare.h>

namespace Pscf {
namespace Rp {

   /**
   * Evaluate the derivative of H with respect to cubic box length L.
   *
   * \see \ref rp_CubicLengthDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class CubicLengthDerivative : public AverageAnalyzer<D,T>
   {
   public:

      /**
      * Constructor.
      */
      CubicLengthDerivative(Simulator<D,T>& simulator, 
                            System<D,T>& system);

      /**
      * Destructor.
      */
      ~CubicLengthDerivative() = default;

   protected:

      /**
      * Compute and return the derivative of H w/ respect to L.
      */
      double compute() override;

      using AverageAnalyzerT = AverageAnalyzer<D,T>;
      using AverageAnalyzer<D,T>::simulator;
      using AverageAnalyzer<D,T>::system;
      using AverageAnalyzer<D,T>::outputFile_;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(CubicLengthDerivative)

} // namespace Rp
} // namespace Pscf
#endif
