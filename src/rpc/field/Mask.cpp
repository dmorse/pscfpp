/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <rpc/field/FieldIo.h>
#include <pscf/cpu/Reduce.h>

#include <rp/field/Mask.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::Mask< 1, Rpc::Types<1> >;
      template class Rp::Mask< 2, Rpc::Types<2> >;
      template class Rp::Mask< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
   }
}
