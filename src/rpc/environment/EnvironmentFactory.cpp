/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EnvironmentFactory.h"  
#include "FilmEnvironment.h"
#include <rp/environment/EnvironmentFactory.tpp>

#if 0
namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   Rp::EnvironmentFactory<D, Rpc::Types<D> >::EnvironmentFactory(Rp::System<D, Rpc::Types<D> >& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Environment subclass className.
   */
   template <int D>
   Environment* 
   Rp::EnvironmentFactory<D, Rpc::Types<D> >::factory(std::string const& className) const
   {
      Environment* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;
 
      // Try to match classname
      if (className == "FilmEnvironment") {
         ptr = new Rp::FilmEnvironment<D, Rpc::Types<D> >(*sysPtr_);
      }

      return ptr;
   }

}
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EnvironmentFactory<1, Rpc::Types<1> >;
      template class EnvironmentFactory<2, Rpc::Types<2> >;
      template class EnvironmentFactory<3, Rpc::Types<3> >;
   }
}
