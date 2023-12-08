#ifndef PSPC_BD_STEP_FACTORY_TPP
#define PSPC_BD_STEP_FACTORY_TPP

#include "BdStepFactory.h"  
#include <pspc/simulate/BdSimulator.h>

// Subclasses of BdStep 
//#include "RealMove.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   BdStepFactory<D>::BdStepFactory(BdSimulator<D>& bdSimulator)
    : bdSimulatorPtr_(&bdSimulator)
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
      if (className == "RealMove") {
         //ptr = new RealMove<D>(*bdSimulatorPtr_);
         ptr = 0;
      }

      return ptr;
   }

}
}
#endif
