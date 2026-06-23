#ifndef RP_CUBIC_LENGTH_DERIVATIVE_H
#define RP_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h> // base class template

namespace Pscf {
namespace Rp {

   /**
   * Evaluate the derivative of H with respect to cubic box length L.
   *
   * \see \ref rp_CubicLengthDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class CubicLengthDerivative : public T::AverageAnalyzer
   {

   protected:

      /**
      * Constructor.
      */
      CubicLengthDerivative(typename T::Simulator& simulator, 
                          System<D,T>& system);

      /**
      * Destructor.
      */
      ~CubicLengthDerivative() = default;

      /**
      * Compute and return the derivative of H w/ respect to L.
      */
      double compute() override;

      using AverageAnalyzerT = typename T::AverageAnalyzer;
      using AverageAnalyzerT::simulator;
      using AverageAnalyzerT::system;
      using AverageAnalyzerT::outputFile_;

   };

}
}
#endif
