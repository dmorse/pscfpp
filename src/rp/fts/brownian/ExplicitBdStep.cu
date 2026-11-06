/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecRandom.h>
#include <pscf/backend/cuda/VecOp.h>

#include <rp/fts/brownian/ExplicitBdStep.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::ExplicitBdStep<1,CUT>;
      template class Rp::ExplicitBdStep<2,CUT>;
      template class Rp::ExplicitBdStep<3,CUT>;
   }
}
