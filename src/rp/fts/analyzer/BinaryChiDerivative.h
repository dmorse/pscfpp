#ifndef RP_CHI_DERIVATIVE_H
#define RP_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h> // base class template

#include <pscf/backends/TmplDeclare.h>
#include <iostream>

namespace Pscf {
namespace Rp {

   /**
   * Evaluate derivative of H with respect to chi for binary system.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_BinaryChiDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class BinaryChiDerivative : public AverageAnalyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryChiDerivative(
                    Simulator<D,T>& simulator, 
                    System<D,T>& system);

      /**
      * Destructor.
      */
      ~BinaryChiDerivative() = default;

   protected:

      /**
      * Compute and return the derivative of H w/ respect to chi.
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
   PSCF_TMPL_DECLARE(BinaryChiDerivative)

} // namespace Rp
} // namespace Pscf
#endif
