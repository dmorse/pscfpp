/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdMove.h"                         // class header
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/VecOp.h>

#include <rp/fts/montecarlo/BdMove.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::BdMove<1, Rpg::Types<1> >;
      template class Rp::BdMove<2, Rpg::Types<2> >;
      template class Rp::BdMove<3, Rpg::Types<3> >;
   }
}
