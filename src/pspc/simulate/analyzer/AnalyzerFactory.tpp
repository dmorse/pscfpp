#ifndef PSPC_MC_MOVE_FACTORY_TPP
#define PSPC_MC_MOVE_FACTORY_TPP

#include "AnalyzerFactory.h"  

// Subclasses of Analyzer 
#include <pspc/simulate/McSimulator.h>
#include "TrajectoryWriter.h"
#include "HamiltonianAnalyzer.h"
#include "BinaryStructureFactorGrid.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AnalyzerFactory<D>::AnalyzerFactory(McSimulator<D>& mcSimulator, System<D>& system)
    : sysPtr_(&system),
      mcSimulatorPtr_(&mcSimulator)
   {}

   /* 
   * Return a pointer to a instance of Analyzer subclass className.
   */
   template <int D>
   Analyzer<D>* AnalyzerFactory<D>::factory(const std::string &className) const
   {
      Analyzer<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      
      // Try to match classname
      if (className == "TrajectoryWriter") {
         ptr = new TrajectoryWriter<D>(*mcSimulatorPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new HamiltonianAnalyzer<D>(*mcSimulatorPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactorGrid") {
         ptr = new BinaryStructureFactorGrid<D>(*mcSimulatorPtr_, *sysPtr_);
      }

      return ptr;
   }

}
}
#endif
