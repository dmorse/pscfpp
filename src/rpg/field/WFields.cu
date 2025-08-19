/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.tpp"

namespace Pscf {
namespace Prdc {
   template class WFieldsTmpl<1, RField<1>, Rpg::FieldIo<1> >;
   template class WFieldsTmpl<2, RField<2>, Rpg::FieldIo<2> >;
   template class WFieldsTmpl<3, RField<3>, Rpg::FieldIo<3> >;
}
namespace Rpg {
   template class WFields<1>;
   template class WFields<2>;
   template class WFields<3>;
} // namespace Rpg
} // namespace Pscf
