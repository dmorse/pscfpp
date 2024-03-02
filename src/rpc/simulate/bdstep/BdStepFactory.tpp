#ifndef RPC_BD_STEP_FACTORY_TPP
#define RPC_BD_STEP_FACTORY_TPP

#include "BdStepFactory.h"  
#include <rpc/simulate/bdstep/BdSimulator.h>

// Subclasses of BdStep 
#include "ExplicitBdStep.h"
#include "MidstepBdStep.h"
#include "PredCorrBdStep.h"
#include "LMBdStep.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   BdStepFactory<D>::BdStepFactory(BdSimulator<D>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of BdStep subclass className.
   */
   template <int D>
   BdStep<D>* BdStepFactory<D>::factory(const std::string &className) const
   {
      BdStep<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      
      // Try to match classname
      if (className == "ExplicitBdStep" || className == "BdStep") {
         ptr = new ExplicitBdStep<D>(*simulatorPtr_);
      } else
      if (className == "MidstepBdStep") {
         ptr = new MidstepBdStep<D>(*simulatorPtr_);
      } else
      if (className == "PredCorrBdStep") {
         ptr = new PredCorrBdStep<D>(*simulatorPtr_);
      } else
      if (className == "LMBdStep") {
         ptr = new LMBdStep<D>(*simulatorPtr_);
      }

      return ptr;
   }

}
}
#endif
