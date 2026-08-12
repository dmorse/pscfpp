/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"                // class header
#include <cp/field/WFields.tpp>     // base class implementation

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   template <int D>
   void WFields<D>::assignField(CField<D,CPT>& lhs, 
                                CField<D,CPT> const & rhs) const
   {
      int n = rhs.capacity();
      UTIL_CHECK(lhs.capacity() == n);
      UTIL_CHECK(meshSize() == n);
      for (int i = 0; i < n; ++i) {
         lhs[i][0] = rhs[i][0];
         lhs[i][1] = rhs[i][1];
      }
   }

} // namespace Cpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Cp {
      template class WFields<1, Prdc::CField<1,CPT>, Cpc::FieldIo<1> >;
      template class WFields<2, Prdc::CField<2,CPT>, Cpc::FieldIo<2> >;
      template class WFields<3, Prdc::CField<3,CPT>, Cpc::FieldIo<3> >;
   }
   namespace Cpc {
      template class WFields<1>;
      template class WFields<2>;
      template class WFields<3>;
   }
}
