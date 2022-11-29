#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include "NrIterator.h"
#include "AmIterator.h"
#include "BinaryRelaxIterator.h"

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
      } else if (className == "BinaryRelaxIterator") {
         ptr = new BinaryRelaxIterator(*sysPtr_);
      } else if (className == "AmIterator") {
         ptr = new AmIterator(*sysPtr_);
      }

      return ptr;
   }

}
}
