/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"               // header

#include <rpg/field/FieldIo.h>
#include <pscf/cuda/Reduce.h>

#include <rp/field/Mask.tpp>    // class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Mask< 1, CudaTp<1> >;
      template class Mask< 2, CudaTp<2> >;
      template class Mask< 3, CudaTp<3> >;
   }
}
