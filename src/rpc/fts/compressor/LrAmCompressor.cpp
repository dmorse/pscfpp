/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <rp/fts/compressor/LrAmCompressor.tpp>

namespace Pscf {
namespace Rpc {

   /*
   * Constructor.
   */
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(Rp::System<D, Rpc::Types<D> >& system)
    : Rp::LrAmCompressor<D, Types<D>, DArray<double> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrAmCompressor<1, Rpc::Types<1>, DArray<double> >;
      template class LrAmCompressor<2, Rpc::Types<2>, DArray<double> >;
      template class LrAmCompressor<3, Rpc::Types<3>, DArray<double> >;
   }
   namespace Rpc {
      template class LrAmCompressor<1>;
      template class LrAmCompressor<2>;
      template class LrAmCompressor<3>;
   }
}
