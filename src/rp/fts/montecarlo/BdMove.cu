/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
#include "BdMove.h"                         // class header
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#endif

#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/backends/CUT.h>

#include <rp/fts/montecarlo/BdMove.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::BdMove<1,CUT>;
      template class Rp::BdMove<2,CUT>;
      template class Rp::BdMove<3,CUT>;
   }
}
