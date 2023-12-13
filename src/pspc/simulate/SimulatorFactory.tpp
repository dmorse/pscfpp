#ifndef PSPC_SIMULATOR_FACTORY_TPP
#define PSPC_SIMULATOR_FACTORY_TPP

#include "SimulatorFactory.h"  
#include <pspc/System.h>

// Subclasses of Simulator 
#include <pspc/simulate/mcmove/McSimulator.h>
#include <pspc/simulate/bdstep/BdSimulator.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   SimulatorFactory<D>::SimulatorFactory(System<D>& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Simulator subclass className.
   */
   template <int D>
   Simulator<D>* SimulatorFactory<D>::factory(const std::string &className) 
   const
   {
      Simulator<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
      
      // Try to match classname
      if (className == "McSimulator" || className == "Simulator") {
         ptr = new McSimulator<D>(*systemPtr_);
      } else
      if (className == "BdSimulator") {
         ptr = new BdSimulator<D>(*systemPtr_);
      } 

      return ptr;
   }

}
}
#endif
