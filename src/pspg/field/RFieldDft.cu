#define PSPG_R_FIELD_DFT_CU

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.tpp"

namespace Pscf {
namespace Pspg {

   template class RFieldDft<1>;
   template class RFieldDft<2>;
   template class RFieldDft<3>;

} // namespace Pscf::Pspg
} // namespace Pscf
