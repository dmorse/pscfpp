
/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"                       // header
#include "McMove.h"
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/CudaVecRandom.h>

#include <rp/fts/montecarlo/RealMove.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RealMove<1, Rpg::Types<1> >;
      template class RealMove<2, Rpg::Types<2> >;
      template class RealMove<3, Rpg::Types<3> >;
   }
}
