#ifndef PSPG_ITERATOR_FACTORY_TPP
#define PSPG_ITERATOR_FACTORY_TPP

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include <pscf/iterator/AmIterator.h>
#include <pspg/iterator/AmStrategyCUDA.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D>::IteratorFactory(IteratorMediatorCUDA<D>& iterMed)
    : iterMedPtr_(&iterMed)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<FieldCUDA>* 
   IteratorFactory<D>::factory(const std::string &className) const
   {
      Iterator<FieldCUDA>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "AmIterator") {
         AmStrategyCUDA* stratPtr= new AmStrategyCUDA();
         ptr = new AmIterator<FieldCUDA>(*iterMedPtr_, *stratPtr);
      }

      return ptr;
   }

}
}
#endif
