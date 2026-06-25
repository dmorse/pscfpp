/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LMBdStep.h"                     // class header
#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/VecOp.h>

#include <rp/fts/brownian/LMBdStep.tpp>   // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   LMBdStep<D>::LMBdStep(Rp::BdSimulator<D, Rpg::Types<D> >& simulator)
    : Rp::LMBdStep<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::LMBdStep<1, Rpg::Types<1> >;
      template class Rp::LMBdStep<2, Rpg::Types<2> >;
      template class Rp::LMBdStep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class LMBdStep<1>;
      template class LMBdStep<2>;
      template class LMBdStep<3>;
   }
}
