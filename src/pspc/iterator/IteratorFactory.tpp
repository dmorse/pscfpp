#ifndef PSPC_ITERATOR_FACTORY_TPP
#define PSPC_ITERATOR_FACTORY_TPP

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include <pspc/iterator/AmIterator.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D>::IteratorFactory(System<D>& system)
    : sys_(&system)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<D>* 
   IteratorFactory<D>::factory(const std::string &className) const
   {
      Iterator<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "AmIterator") {
         ptr = new AmIterator<D>(*sys_);
      }

      return ptr;
   }

}
}
#endif
