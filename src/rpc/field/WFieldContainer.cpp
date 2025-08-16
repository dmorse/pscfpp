/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldContainer.tpp"

namespace Pscf {
namespace Prdc {
   template class WFieldsTmpl<1, RField<1>, Rpc::FieldIo<1> >;
   template class WFieldsTmpl<2, RField<2>, Rpc::FieldIo<2> >;
   template class WFieldsTmpl<3, RField<3>, Rpc::FieldIo<3> >;
}
namespace Rpc {

   template class WFieldContainer<1>;
   template class WFieldContainer<2>;
   template class WFieldContainer<3>;

} // namespace Rpc
} // namespace Pscf
