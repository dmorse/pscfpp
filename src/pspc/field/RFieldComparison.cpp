#define PSPC_R_FIELD_COMPARISON_CPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.h"

namespace Pscf {
namespace Pspc
{

   template class RFieldComparison<1>;
   template class RFieldComparison<2>;
   template class RFieldComparison<3>;

} // namespace Pscf::Pspc
} // namespace Pscf
