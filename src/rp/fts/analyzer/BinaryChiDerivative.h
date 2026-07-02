#ifndef RP_CHI_DERIVATIVE_H
#define RP_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h> // base class template

#include <iostream>

namespace Pscf {
namespace Rp {

   /**
   * Evaluate derivative of H with respect to chi for binary system.
   *
   * Specializations of this template are used as base classes for two
   * closely analogous class templates, also named BinaryChiDerivative, that 
   * are defined in the Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>.
   *
   * \see \ref rp_BinaryChiDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class BinaryChiDerivative : public AverageAnalyzer<D,T>
   {

   protected:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryChiDerivative(Simulator<D,T>& simulator, 
                    System<D,T>& system);

      /**
      * Destructor.
      */
      ~BinaryChiDerivative() = default;

      /**
      * Compute and return the derivative of H w/ respect to chi.
      */
      double compute() override;

      // Aliases for base classes.
      using AnalyzerT = Analyzer<D,T>;
      using AverageAnalyzerT = AverageAnalyzer<D,T>;

      // Inherited protected member functions (selected).
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   };

}
}
#endif
