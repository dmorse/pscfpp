#ifndef RPC_W_FIELD_CONTAINER_TPP
#define RPC_W_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldContainer.h"           // class header
#include <prdc/field/WFieldsReal.tpp>  // base class implementation

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   template <int D>
   void WFieldContainer<D>::assignRField(RField<D>& lhs, 
                                         RField<D> const & rhs) const
   {
      int n = rhs.capacity();
      UTIL_CHECK(lhs.capacity() == n);
      UTIL_CHECK(meshSize() == n);
      for (int i = 0; i < n; ++i) {
         lhs[i] = rhs[i];
      }
   }

} // namespace Rpc
} // namespace Pscf
#endif
