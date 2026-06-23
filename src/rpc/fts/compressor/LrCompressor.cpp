/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "LrCompressor.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <prdc/cpu/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>

#include <rp/fts/compressor/LrCompressor.tpp>

namespace Pscf {
namespace Rpc {

   /*
   * Constructor.
   */
   template <int D>
   LrCompressor<D>::LrCompressor(Rp::System<D, Rpc::Types<D> >& system)
    : Rp::LrCompressor<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrCompressor<1, Rpc::Types<1> >;
      template class LrCompressor<2, Rpc::Types<2> >;
      template class LrCompressor<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class LrCompressor<1>;
      template class LrCompressor<2>;
      template class LrCompressor<3>;
   }
}
