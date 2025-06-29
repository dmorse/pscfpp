/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldContainer.tpp"

namespace Pscf {
namespace Prdc {
   template class WFieldsReal<1, RField<1>, Rpg::FieldIo<1> >;
   template class WFieldsReal<2, RField<2>, Rpg::FieldIo<2> >;
   template class WFieldsReal<3, RField<3>, Rpg::FieldIo<3> >;
}
namespace Rpg {
   template class WFieldContainer<1>;
   template class WFieldContainer<2>;
   template class WFieldContainer<3>;
} // namespace Rpg
} // namespace Pscf
