/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CpuVecRandom.h>
#include <pscf/backend/cpp/VecOp.h>

#include <rp/fts/brownian/LMBdStep.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LMBdStep<1,CPT>;
      template class LMBdStep<2,CPT>;
      template class LMBdStep<3,CPT>;
   }
}
