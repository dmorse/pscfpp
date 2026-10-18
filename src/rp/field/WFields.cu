/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/DeviceArray.h>

#include <rp/field/WFields.tpp>   
//#include <rp/field/WFields_u.h>       // class specialziation

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Explicit instantiation definitions - this class
   template class WFields<1,CUT>;
   template class WFields<2,CUT>;
   template class WFields<3,CUT>;

} // namespace Rp
} // namespace Pscf

