/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/brownian/PredCorrBdStep.tpp>  // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PredCorrBdStep<1,CPT>;
      template class PredCorrBdStep<2,CPT>;
      template class PredCorrBdStep<3,CPT>;
   }
}
