#ifndef RPG_MC_MOVE_FACTORY_TPP
#define RPG_MC_MOVE_FACTORY_TPP

#include "AnalyzerFactory.h"  

// Subclasses of Analyzer 
#include "TrajectoryWriter.h"
#include "HamiltonianAnalyzer.h"
#include "BinaryStructureFactorGrid.h"
#include "StepLogger.h"
#include "PerturbationDerivative.h"
#include "ChiDerivative.h"
#include "ConcentrationDerivative.h"
#include "MaxOrderParameter.h"
#include "FourthOrderParameter.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AnalyzerFactory<D>::AnalyzerFactory(Simulator<D>& simulator, System<D>& system)
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
         ptr = 
           new HamiltonianAnalyzer<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactorGrid") {
         ptr = 
           new BinaryStructureFactorGrid<D>(*simulatorPtr_, 
                                            *sysPtr_);
      } else if (className == "StepLogger") {
         ptr = new StepLogger<D>();
      } else if (className == "PerturbationDerivative") {
         ptr = new PerturbationDerivative<D>(*simulatorPtr_, 
                                             *sysPtr_);
      } else if (className == "ChiDerivative") {
         ptr = new ChiDerivative<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "ConcentrationDerivative") {
         ptr = new ConcentrationDerivative<D>(*simulatorPtr_, 
                                              *sysPtr_);
      } else if (className == "MaxOrderParameter") {
         ptr = new MaxOrderParameter<D>(*simulatorPtr_, *sysPtr_);
      } else if (className == "FourthOrderParameter") {
         ptr = new FourthOrderParameter<D>(*simulatorPtr_, *sysPtr_);
      }

      return ptr;
   }

}
}
#endif
