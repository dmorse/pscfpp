#define PSPC_R_FIELD_COMPARISON_CPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.h"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   template class RFieldComparison<1>;
   template class RFieldComparison<2>;
   template class RFieldComparison<3>;

} // namespace Pscf::Prdc::Gpu
} // namespace Pscf::Prdc
} // namespace Pscf
