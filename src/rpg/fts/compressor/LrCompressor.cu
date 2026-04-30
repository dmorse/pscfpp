/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrCompressor.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/fts/compressor/LrCompressor.tpp>

namespace Pscf {
namespace Rpg {

   /*
   * Constructor.
   */
   template <int D>
   LrCompressor<D>::LrCompressor(System<D>& system)
    : Rp::LrCompressor<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrCompressor<1, Rpg::Types<1> >;
      template class LrCompressor<2, Rpg::Types<2> >;
      template class LrCompressor<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class LrCompressor<1>;
      template class LrCompressor<2>;
      template class LrCompressor<3>;
   }
}
