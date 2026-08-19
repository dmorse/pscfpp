/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IteratorFactory.h"

// Subclasses of Iterator
#include <rp/scft/iterator/AmIteratorBasis.h>
#include <rp/scft/iterator/AmIteratorGrid.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   IteratorFactory<D,T>::IteratorFactory(System<D,T>& system)
    : sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D, class T>
   Iterator<D,T>*
   IteratorFactory<D,T>::factory(const std::string &className) const
   {
      Iterator<D,T>* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "Iterator" || className == "AmIteratorBasis"
          || className == "AmIterator" ) {
         ptr = new AmIteratorBasis<D,T>(*sysPtr_);
      } else
      if (className == "AmIteratorGrid") {
         ptr = new AmIteratorGrid<D,T>(*sysPtr_);
      }

      return ptr;
   }

}
}
