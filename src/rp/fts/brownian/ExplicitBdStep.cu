/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/VecOp.h>

#include <rp/fts/brownian/ExplicitBdStep.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::ExplicitBdStep<1,CUT>;
      template class Rp::ExplicitBdStep<2,CUT>;
      template class Rp::ExplicitBdStep<3,CUT>;
   }
}
