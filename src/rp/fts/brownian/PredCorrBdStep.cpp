/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/PredCorrBdStep.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PredCorrBdStep<1,CPT>;
      template class PredCorrBdStep<2,CPT>;
      template class PredCorrBdStep<3,CPT>;
   }
}
