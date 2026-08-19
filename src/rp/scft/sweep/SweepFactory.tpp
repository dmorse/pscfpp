/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"

// Subclasses of Sweep
#include "LinearSweep.h"

namespace Pscf {
namespace Rp {

   using namespace Util;

   template <int D, class T>
   SweepFactory<D,T>::SweepFactory(System<D,T>& system)
    : systemPtr_(&system)
   {}

   /**
   * Return a pointer to a Sweep subclass with name className
   */
   template <int D, class T>
   Sweep<D,T>* 
   SweepFactory<D,T>::factory(std::string const & className) const
   {
       Sweep<D,T>* ptr = nullptr;

       // First check if name is known by any subfactories
       ptr = trySubfactories(className);
       if (ptr) return ptr;

       // Explicit class names
       if (className == "Sweep" || className == "LinearSweep") {
           ptr = new LinearSweep<D,T>(*systemPtr_);
       }

       return ptr;
   }

} // namespace Rp
} // namespace Pscf

