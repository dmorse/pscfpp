/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <prdc/cpu/RField.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/fts/compressor/AmCompressor.tpp>

namespace Pscf {
namespace Rpc {

   /*
   * Constructor.
   */
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
    : Rp::AmCompressor<D, Types<D>, DArray<double> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AmCompressor<1, Rpc::Types<1>, DArray<double> >;
      template class AmCompressor<2, Rpc::Types<2>, DArray<double> >;
      template class AmCompressor<3, Rpc::Types<3>, DArray<double> >;
   }
   namespace Rpc {
      template class AmCompressor<1>;
      template class AmCompressor<2>;
      template class AmCompressor<3>;
   }
}
