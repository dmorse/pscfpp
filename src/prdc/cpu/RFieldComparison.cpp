/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
