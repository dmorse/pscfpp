/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"  

// Subclasses of Sweep 
// #include "CompositionSweep.h"

namespace Pscf {
namespace Fd1d {

   using namespace Util;

   /* 
   * Return a pointer to a instance of Sweep subclass speciesName.
   */
   Sweep* SweepFactory::factory(const std::string &className) const
   {
      Sweep *ptr = 0;

      ptr = trySubfactories(className);
      if (ptr) return ptr;     

      // Try explicit class names
      if (className == "CompositionSweep") {
         ptr = 0;
         // ptr = new Sweep();
      } 

      return ptr;
   }

}
}
