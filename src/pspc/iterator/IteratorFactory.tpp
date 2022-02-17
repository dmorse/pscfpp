#ifndef PSPC_ITERATOR_FACTORY_TPP
#define PSPC_ITERATOR_FACTORY_TPP

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include <pscf/iterator/AmIterator.h>
#include <pspc/iterator/AmStrategyCPU.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D>::IteratorFactory(IteratorMediatorCPU<D>& iterMed)
    : iterMedPtr_(&iterMed)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<FieldCPU>* 
   IteratorFactory<D>::factory(const std::string &className) const
   {
      Iterator<FieldCPU>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "AmIterator") {
         AmStrategyCPU* stratPtr= new AmStrategyCPU();
         ptr = new AmIterator<FieldCPU>(*iterMedPtr_, *stratPtr);
      }

      return ptr;
   }

}
}
#endif
