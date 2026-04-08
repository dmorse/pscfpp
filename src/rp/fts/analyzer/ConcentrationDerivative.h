#ifndef RP_CONCENTRATION_DERIVATIVE_H
#define RP_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Rp {

   /**
   * Evaluate the derivative of H with respect to concentration.
   *
   * Specializations of this template are used as base classes for two
   * closely analogous class templates, also named ConcentrationDerivative, 
   * that are defined in the Rpc and Rpg namespaces for use in the 
   * pscf_rpc and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>.
   *
   * \see \ref rp_ConcentrationDerivative_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class ConcentrationDerivative : public T::AverageAnalyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationDerivative(typename T::Simulator& simulator,
                              typename T::System& system);

   protected:

      /**
      * Compute and return the derivative of H w/ respect to concentration.
      */
      double compute() override;

      #if 0
      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      void outputValue(int step, double value) override;
      #endif

      using AnalyzerT = typename T::Analyzer;
      using AverageAnalyzerT = typename T::AverageAnalyzer;
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   };

} // namespace Rp
} // namespace Pscf
#endif
