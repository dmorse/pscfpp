/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"

#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/math/IntVec.h>

#include <rp/fts/brownian/PredCorrBdStep.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PredCorrBdStep<1,CUT>;
      template class PredCorrBdStep<2,CUT>;
      template class PredCorrBdStep<3,CUT>;
   }
}
