#ifndef CPC_STEP_FACTORY_TPP
#define CPC_STEP_FACTORY_TPP

#include "StepFactory.h"  
#include <cpc/fts/Simulator.h>

// Subclasses of Step 
//#include "ExplicitStep.h"

namespace Pscf {
namespace Cpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   StepFactory<D>::StepFactory(Simulator<D>& simulator)
    : simulatorPtr_(&simulator)
   {}

   /* 
   * Return a pointer to a instance of Step subclass className.
   */
   template <int D>
   Step<D>* StepFactory<D>::factory(const std::string &className) const
   {
      Step<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

     
      #if 0 
      // Try to match classname
      if (className == "ExplicitStep" || className == "Step") {
         ptr = new ExplicitStep<D>(*simulatorPtr_);
      } 
      #endif

      return ptr;
   }

} // namespace Cpc
} // namespace Pscf
#endif
