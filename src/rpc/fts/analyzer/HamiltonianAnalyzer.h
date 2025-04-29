#ifndef RPC_HAMILTONIAN_ANALYZER_H
#define RPC_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Compute averages and output block averages of Hamiltonian components.
   *
   * This class computes separate averages for each component of the
   * total simulation Hamiltonian (ideal gas contributions (lnQ) and 
   * Field contribution (HW)) as well as for the total, and 
   * periodically outputs block averages of each to a file.
   *
   * \see \ref rpc_HamiltonianAnalyzer_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class HamiltonianAnalyzer : public AverageListAnalyzer<D>
   {

   public:
   
      /**
      * Constructor.
      */
      HamiltonianAnalyzer(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~HamiltonianAnalyzer()
      {} 
      
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

   protected:

      using ParamComposite::setClassName;
      using AverageListAnalyzer<D>::setName;
      using AverageListAnalyzer<D>::setValue;
      using AverageListAnalyzer<D>::simulator;
      using AverageListAnalyzer<D>::system;
      
      /**
      * Compute and store values of Hamiltonian components.
      */  
      void compute();
     
   private: 

      /// Array index for ideal gas contributions (lnQ) accumulator.
      int idealId_;

      /// Array index for field contribution (HW) accumulator.
      int fieldId_;

      /// Array index for total Hamiltonian accumulator.
      int totalId_;

   };
  
   #ifndef RPC_HAMILTONIAN_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class HamiltonianAnalyzer<1>;
   extern template class HamiltonianAnalyzer<2>;
   extern template class HamiltonianAnalyzer<3>;
   #endif

}
}
#endif 
