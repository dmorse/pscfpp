/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"  

// Subclasses of Sweep 
#include "LinearSweep.h"

namespace Pscf {
namespace R1d {

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
