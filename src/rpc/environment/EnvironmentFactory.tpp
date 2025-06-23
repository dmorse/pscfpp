#ifndef RPC_ENVIRONMENT_FACTORY_TPP
#define RPC_ENVIRONMENT_FACTORY_TPP

#include "EnvironmentFactory.h"  
#include "MixAndMatchEnvs.h" // Subclasses of Environment

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   EnvironmentFactory<D>::EnvironmentFactory(System<D>& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Environment subclass className.
   */
   template <int D>
   Environment* EnvironmentFactory<D>::factory(const std::string &className) 
   const
   {
      Environment* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "FilmEnvironment") {
         ptr = new FilmEnvironment<D>(*sysPtr_);
      }

      return ptr;
   }

}
}
#endif
