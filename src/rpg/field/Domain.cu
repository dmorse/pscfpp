/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"               // class header
#include <rp/field/FieldIo.h>
#include <prdc/field/cuda/WaveList.h>
#include <prdc/field/cuda/FFT.h>

#include <rp/field/Domain.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Domain<1,CUT>;
      template class Domain<2,CUT>;
      template class Domain<3,CUT>;
   }
}
