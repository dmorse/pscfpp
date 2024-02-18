#ifndef PSPC_ANALYZER_FACTORY_TPP
#define PSPC_ANALYZER_FACTORY_TPP

#include "AnalyzerFactory.h"  

// Subclasses of Analyzer 
#include "TrajectoryWriter.h"
#include "HamiltonianAnalyzer.h"
#include "HamiltonianAutoCorr.h"
#include "BinaryStructureFactorGrid.h"
#include "StepLogger.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AnalyzerFactory<D>::AnalyzerFactory(Simulator<D>& simulator, 
                                       System<D>& system)
    : sysPtr_(&system),
      simulatorPtr_(&simulator)
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
         ptr = new TrajectoryWriter<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new HamiltonianAnalyzer<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "HamiltonianAutoCorr") {
         ptr = new HamiltonianAutoCorr<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactorGrid") {
         ptr 
           = new BinaryStructureFactorGrid<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "StepLogger") {
         ptr = new StepLogger<D>();
      }

      return ptr;
   }

}
}
#endif
