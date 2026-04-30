/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"
#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/math/IntVec.h>

#include <rp/fts/brownian/PredCorrBdStep.tpp>  // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   PredCorrBdStep<D>::PredCorrBdStep(BdSimulator<D>& simulator)
    : Rp::PredCorrBdStep<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::PredCorrBdStep<1, Rpg::Types<1> >;
      template class Rp::PredCorrBdStep<2, Rpg::Types<2> >;
      template class Rp::PredCorrBdStep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class PredCorrBdStep<1>;
      template class PredCorrBdStep<2>;
      template class PredCorrBdStep<3>;
   }
}
