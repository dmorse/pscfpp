/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <pscf/cpu/Reduce.h>

#include <rp/field/Mask.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cpu;
      template class Rp::Mask< 1, RField<1>, Rpc::FieldIo<1> >;
      template class Rp::Mask< 2, RField<2>, Rpc::FieldIo<2> >;
      template class Rp::Mask< 3, RField<3>, Rpc::FieldIo<3> >;
   }
   namespace Rpc {
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
   }
}
