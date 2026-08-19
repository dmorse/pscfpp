/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EnvironmentFactory.h"  
#include "FilmEnvironment.h"

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor
   */
   template <int D, class T>
   EnvironmentFactory<D,T>::EnvironmentFactory(System<D,T>& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Environment subclass className.
   */
   template <int D, class T>
   Environment* 
   EnvironmentFactory<D,T>::factory(std::string const& className) const
   {
      Environment* ptr = nullptr;

      // Try subfactories first
      ptr = Factory<Environment>::trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "FilmEnvironment") {
         ptr = new FilmEnvironment<D,T>(*sysPtr_);
      }

      return ptr;
   }

}
}

