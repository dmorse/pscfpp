/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <pscf/cuda/Reduce.h>

#include <rp/field/Mask.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      template class Mask< 1, RField<1>, Rpg::FieldIo<1> >;
      template class Mask< 2, RField<2>, Rpg::FieldIo<2> >;
      template class Mask< 3, RField<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
   }
}
