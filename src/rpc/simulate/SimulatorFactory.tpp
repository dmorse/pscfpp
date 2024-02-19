#ifndef RPC_SIMULATOR_FACTORY_TPP
#define RPC_SIMULATOR_FACTORY_TPP

#include "SimulatorFactory.h"  
#include <rpc/System.h>

// Subclasses of Simulator 
#include <rpc/simulate/mcmove/McSimulator.h>
#include <rpc/simulate/bdstep/BdSimulator.h>

namespace Pscf {
namespace Rpc {

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
