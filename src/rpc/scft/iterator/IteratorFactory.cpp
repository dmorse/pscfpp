/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IteratorFactory.h"  

// Subclasses of Rpc::Iterator 
#include "AmIteratorBasis.h"
#include "AmIteratorGrid.h"

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   IteratorFactory<D, Rpc::Types<D> >::IteratorFactory(Rp::System<D, Rpc::Types<D> >& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Iterator subclass className.
   */
   template <int D>
   Iterator<D, Rpc::Types<D> >* 
   IteratorFactory<D, Rpc::Types<D> >::factory(const std::string &className) 
   const
   {
      Rp::Iterator<D, Rpc::Types<D> >* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "Iterator" || className == "AmIteratorBasis" 
          || className == "AmIterator" ) {
         ptr = new Rp::AmIteratorBasis<D, Rpc::Types<D> >(*sysPtr_);
      } else 
      if (className == "AmIteratorGrid") {
         ptr = new Rp::AmIteratorGrid<D, Rpc::Types<D> >(*sysPtr_);
      }

      return ptr;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
namespace Rp {
   template class IteratorFactory<1, Rpc::Types<1> >;
   template class IteratorFactory<2, Rpc::Types<2> >;
   template class IteratorFactory<3, Rpc::Types<3> >;
}
}
