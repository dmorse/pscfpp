/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <rpg/field/FieldIo.h>
#include <pscf/cuda/Reduce.h>

#include <rp/field/Mask.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Mask< 1, Rpg::Types<1> >;
      template class Mask< 2, Rpg::Types<2> >;
      template class Mask< 3, Rpg::Types<3> >;
   }
   #if 0
   namespace Rpg {
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
   }
   #endif
}
