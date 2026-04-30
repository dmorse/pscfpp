/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/VecOpMisc.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/math/IntVec.h>

#include <util/global.h>

#include <rp/fts/compressor/LrAmCompressor.tpp>

namespace Pscf {
namespace Rpg {

   /*
   * Constructor.
   */
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
    : Rp::LrAmCompressor<D, Types<D>, DeviceArray<cudaReal> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template 
      class LrAmCompressor<1, Rpg::Types<1>, DeviceArray<cudaReal> >;
      template 
      class LrAmCompressor<2, Rpg::Types<2>, DeviceArray<cudaReal> >;
      template 
      class LrAmCompressor<3, Rpg::Types<3>, DeviceArray<cudaReal> >;
   }
   namespace Rpg {
      template class LrAmCompressor<1>;
      template class LrAmCompressor<2>;
      template class LrAmCompressor<3>;
   }
}
