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
   using namespace Pscf::Prdc::Cpu;

   template <int D>
   void WFields<D>::assignField(CField<D, CppTp<D> >& lhs, 
                                CField<D, CppTp<D> > const & rhs) const
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
      template class WFields<1, Prdc::CField<1, CppTp<1> >, Cpc::FieldIo<1> >;
      template class WFields<2, Prdc::CField<2, CppTp<2> >, Cpc::FieldIo<2> >;
      template class WFields<3, Prdc::CField<3, CppTp<3> >, Cpc::FieldIo<3> >;
   }
   namespace Cpc {
      template class WFields<1>;
      template class WFields<2>;
      template class WFields<3>;
   }
}
