/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"               // template header

#include <rpc/field/FieldIo.h>
#include <pscf/cpu/Reduce.h>

#include <rp/field/Mask.tpp>    // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::Mask< 1, CppTp<1> >;
      template class Rp::Mask< 2, CppTp<2> >;
      template class Rp::Mask< 3, CppTp<3> >;
   }
}
