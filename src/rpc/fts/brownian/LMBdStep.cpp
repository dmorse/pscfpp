/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LMBdStep.h"                     // class header
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/cpu/VecOp.h>

#include <rp/fts/brownian/LMBdStep.tpp>   // base class implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   LMBdStep<D>::LMBdStep(Rp::BdSimulator<D, Rpc::Types<D> >& simulator)
    : Rp::LMBdStep<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::LMBdStep<1, Rpc::Types<1> >;
      template class Rp::LMBdStep<2, Rpc::Types<2> >;
      template class Rp::LMBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class LMBdStep<1>;
      template class LMBdStep<2>;
      template class LMBdStep<3>;
   }
}
