/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include "AmIteratorBasis.h"
#include "AmIteratorGrid.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D>::IteratorFactory(System<D>& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<D>* IteratorFactory<D>::factory(std::string const& className) 
   const
   {
      Iterator<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "Iterator" || className == "AmIteratorBasis") {
         ptr = new AmIteratorBasis<D>(*sysPtr_);
      } else 
      if (className == "AmIteratorGrid") {
         ptr = new AmIteratorGrid<D>(*sysPtr_);
      }

      return ptr;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rpg {
      template class IteratorFactory<1>;
      template class IteratorFactory<2>;
      template class IteratorFactory<3>;
   }
}
