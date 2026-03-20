/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"               // class header
#include <rpg/field/FieldIo.h>
#include <prdc/cuda/WaveList.h>
#include <prdc/cuda/FFT.h>

#include <rp/field/Domain.tpp>    // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      template class Domain<1, FFT<1>, WaveList<1>, Rpg::FieldIo<1> >;
      template class Domain<2, FFT<2>, WaveList<2>, Rpg::FieldIo<2> >;
      template class Domain<3, FFT<3>, WaveList<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      template class Domain<1>;
      template class Domain<2>;
      template class Domain<3>;
   }
}
