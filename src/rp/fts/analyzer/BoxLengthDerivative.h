#ifndef RP_BOX_LENGTH_DERIVATIVE_H
#define RP_BOX_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to cubic box length L.
   *
   * \see \ref rp_BoxLengthDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class BoxLengthDerivative : public T::AverageAnalyzer
   {

   public:

      /**
      * Constructor.
      */
      BoxLengthDerivative(typename T::Simulator& simulator, 
                          typename T::System& system);

   protected:

      using AverageAnalyzerT = typename T::AverageAnalyzer;
      using AverageAnalyzerT::simulator;
      using AverageAnalyzerT::system;
      using AverageAnalyzerT::outputFile_;

      /**
      * Compute and return the derivative of H w/ respect to L.
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

}
}
#endif
