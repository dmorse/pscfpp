/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"               // class header
#include <rpg/field/FieldIo.h>
#include <prdc/field/cuda/WaveList.h>
#include <prdc/field/cuda/FFT.h>

#include <rp/field/Domain.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Domain<1, CudaTp<1> >;
      template class Domain<2, CudaTp<2> >;
      template class Domain<3, CudaTp<3> >;
   }
}
