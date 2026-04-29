/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/math/IntVec.h>

#include <rp/fts/brownian/PredCorrBdStep.tpp>  // base class implementation

namespace Pscf {
namespace Rpc {

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
      template class Rp::PredCorrBdStep<1, Rpc::Types<1> >;
      template class Rp::PredCorrBdStep<2, Rpc::Types<2> >;
      template class Rp::PredCorrBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class PredCorrBdStep<1>;
      template class PredCorrBdStep<2>;
      template class PredCorrBdStep<3>;
   }
}
