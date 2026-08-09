/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdMove.h"                         // class header
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/cpu/VecOp.h>

#include <rp/fts/montecarlo/BdMove.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::BdMove<1,CPT>;
      template class Rp::BdMove<2,CPT>;
      template class Rp::BdMove<3,CPT>;
   }
}
