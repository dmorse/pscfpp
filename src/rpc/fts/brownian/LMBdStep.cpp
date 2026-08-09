/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LMBdStep.h"                     // class header

#include <rp/fts/brownian/BdSimulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/cpu/VecOp.h>

#include <rp/fts/brownian/LMBdStep.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LMBdStep<1,CPT>;
      template class LMBdStep<2,CPT>;
      template class LMBdStep<3,CPT>;
   }
}
