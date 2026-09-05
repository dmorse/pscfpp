/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/CudaVecRandom.h>

#include <rp/fts/brownian/PredCorrBdStep.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PredCorrBdStep<1,CUT>;
      template class PredCorrBdStep<2,CUT>;
      template class PredCorrBdStep<3,CUT>;
   }
}
