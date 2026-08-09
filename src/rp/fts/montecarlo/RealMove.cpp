
/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
#include "RealMove.h"                       // header
#include "McMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#endif

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/backends/CPT.h>

#include <rp/fts/montecarlo/RealMove.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RealMove<1,CPT>;
      template class RealMove<2,CPT>;
      template class RealMove<3,CPT>;
   }
}
