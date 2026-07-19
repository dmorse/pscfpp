
/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"                       // header
#include "McMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/montecarlo/RealMove.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RealMove<1, Cpp<1> >;
      template class RealMove<2, Cpp<2> >;
      template class RealMove<3, Cpp<3> >;
   }
}
