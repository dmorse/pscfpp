/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecRandom.h>
#include <pscf/backend/cuda/VecOp.h>

#include <rp/fts/montecarlo/BdMove.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdMove<1,CUT>;
      template class BdMove<2,CUT>;
      template class BdMove<3,CUT>;
   }
}
