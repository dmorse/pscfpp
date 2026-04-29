#ifndef RP_HAMILTONIAN_ANALYZER_H
#define RP_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace Pscf {
namespace Rp {

   /**
   * Compute averages and output block averages of Hamiltonian components.
   *
   * This class computes separate averages for each component of the
   * total simulation Hamiltonian (ideal gas contributions (lnQ) and
   * Field contribution (HW)) as well as for the total, and
   * periodically outputs block averages of each to a file.
   *
   * Specializations of this template are used as base classes for two
   * closely analogous class templates, also named HamiltonianAnalyzer, 
   * that are defined in the Rpc and Rpg namespaces for use in the 
   * pscf_rpc and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_HamiltonianAnalyzer_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class HamiltonianAnalyzer : public T::AverageListAnalyzer
   {

   public:

      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      void readParameters(std::istream& in) override;

   protected:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      HamiltonianAnalyzer(typename T::Simulator& simulator, 
                          typename T::System& system);

      /**
      * Destructor.
      */
      ~HamiltonianAnalyzer() = default;
     
      /**
      * Compute and store values of Hamiltonian components.
      */
      void compute() override;

      // Aliases for base classes.
      using AverageListAnalyzerT = typename T::AverageListAnalyzer;
      using AnalyzerT = typename T::Analyzer;

      // Inherited protected member functions (selected).
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   private:

      /// Array index for ideal gas contributions (lnQ) accumulator.
      int idealId_;

      /// Array index for field contribution (HW) accumulator.
      int fieldId_;

      /// Array index for total Hamiltonian accumulator.
      int totalId_;

   };

}
}
#endif
