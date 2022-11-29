/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"  

// Subclasses of Sweep 
#include "LinearSweep.h"

namespace Pscf {
namespace Fd1d {

   using namespace Util;

   SweepFactory::SweepFactory(System& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Sweep subclass speciesName.
   */
   Sweep* SweepFactory::factory(const std::string &className) const
   {
      Sweep *ptr = 0;

      // First if name is known by any subfactories
      ptr = trySubfactories(className);
      if (ptr) return ptr;     

      // Explicit class names
      if (className == "Sweep" || className == "LinearSweep") {
         ptr = new LinearSweep(*systemPtr_);
      } 

      return ptr;
   }

}
}
