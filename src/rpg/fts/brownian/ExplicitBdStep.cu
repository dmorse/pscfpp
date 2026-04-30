/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExplicitBdStep.h"                   // class header
#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <rpg/field/CFields.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/VecOp.h>

#include <rp/fts/brownian/ExplicitBdStep.tpp> // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   ExplicitBdStep<D>::ExplicitBdStep(BdSimulator<D>& simulator)
    : Rp::ExplicitBdStep<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::ExplicitBdStep<1, Rpg::Types<1> >;
      template class Rp::ExplicitBdStep<2, Rpg::Types<2> >;
      template class Rp::ExplicitBdStep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class ExplicitBdStep<1>;
      template class ExplicitBdStep<2>;
      template class ExplicitBdStep<3>;
   }
}
