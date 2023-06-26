#ifndef PSPC_MC_HAMILTONIAN_ANALYZER_H
#define PSPC_MC_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"
#include <pspc/System.h>

namespace Util{
   class Average;
}

namespace Pscf {
namespace Pspc 
{
   using namespace Util;

   /**
   * Compute averages and output block averages of Hamiltonian components.
   *
   * This class computes separate averages for each component of the
   * total simulation Hamiltonian (ideal gas contributions (lnQ) and 
   * Field contribution (HW)) as well as for the total, and 
   * periodically outputs block averages of each to a file.
   *
   * \ingroup Pspc_Analyzer_Module
   */
   template <int D>
   class McHamiltonianAnalyzer : public AverageListAnalyzer<D>
   {

   public:
   
      /**
      * Constructor.
      */
      McHamiltonianAnalyzer(McSimulator<D>& mcSimulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~McHamiltonianAnalyzer()
      {} 
      
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

      #if 0
      /**
      * Load parameters from archive when restarting. 
      *
      * \param ar loading/input archive
      */
      virtual void loadParameters(Serializable::IArchive& ar); 
   
      /**
      * Save internal state to archive.
      *
      * \param ar saving/output archive
      */
      virtual void save(Serializable::OArchive& ar);
      #endif

      /**
      * Pointer to parent Simulator
      */
      McSimulator<D>* mcSimulatorPtr_;     
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_; 
      
      using ParamComposite::setClassName;
      using AverageListAnalyzer<D>::setName;
      using AverageListAnalyzer<D>::setValue;
      
   protected:

      /**
      * Compute and store values of Hamiltonian components.
      */  
      void compute();
      
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /** 
      * Return reference to parent McSimulator.
      */
      McSimulator<D>& mcSimulator();
 
   private: 
 
      /// Array index for ideal gas contributions (lnQ) accumulator.
      int idealId_;

      /// Array index for field contribution (HW) accumulator.
      int fieldId_;

      /// Array index for total Hamiltonian accumulator.
      int totalId_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& McHamiltonianAnalyzer<D>::system()
   {  return *systemPtr_; }
   
   //Get parent McSimulator object.
   template <int D>
   inline McSimulator<D>& McHamiltonianAnalyzer<D>::mcSimulator()
   {  return *mcSimulatorPtr_; }


}
}
#endif 
