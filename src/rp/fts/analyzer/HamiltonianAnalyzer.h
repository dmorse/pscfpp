#ifndef RP_HAMILTONIAN_ANALYZER_H
#define RP_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageListAnalyzer.h> // base class template

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
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Cpp<D> or CudaTp<D>
   *
   * \see \ref rp_HamiltonianAnalyzer_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class HamiltonianAnalyzer : public AverageListAnalyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      HamiltonianAnalyzer(Simulator<D,T>& simulator, 
                          System<D,T>& system);

      /**
      * Destructor.
      */
      ~HamiltonianAnalyzer() = default;
     
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      void readParameters(std::istream& in) override;

   protected:

      /**
      * Compute and store values of Hamiltonian components.
      */
      void compute() override;

      // Aliases for base classes.
      using AverageListAnalyzerT = AverageListAnalyzer<D,T>;
      using AnalyzerT = Analyzer<D,T>;

      // Inherited protected member functions (selected).
      using Analyzer<D,T>::simulator;
      using Analyzer<D,T>::system;

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
