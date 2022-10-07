#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include "NrIterator.h"
#include "FdIterator.h"
#include "AmIterator.h"

namespace Pscf {
namespace Fd1d {

   using namespace Util;

   /*
   * Constructor
   */
   IteratorFactory::IteratorFactory(System& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   Iterator* IteratorFactory::factory(const std::string &className) const
   {
      Iterator* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "Iterator" || className == "NrIterator") {
         ptr = new NrIterator(*sysPtr_);
      } else if (className == "FdIterator") {
         ptr = new FdIterator(*sysPtr_);
      } else if (className == "AmIterator") {
         ptr = new AmIterator(*sysPtr_);
      }

      return ptr;
   }

}
}
