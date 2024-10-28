#ifndef RPC_ITERATOR_FACTORY_TPP
#define RPC_ITERATOR_FACTORY_TPP

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include "AmIteratorBasis.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D>::IteratorFactory(System<D>& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<D>* IteratorFactory<D>::factory(const std::string &className) const
   {
      Iterator<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "Iterator" || className == "AmIteratorBasis" 
                                  || className == "AmIterator" ) {
         ptr = new AmIteratorBasis<D>(*sysPtr_);
      }

      return ptr;
   }

}
}
#endif
