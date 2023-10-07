#define PSPG_R_FIELD_CU

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.tpp"

namespace Pscf {
namespace Pspg {

   template class RField<1>;
   template class RField<2>;
   template class RField<3>;

} // namespace Pscf::Pspg
} // namespace Pscf
