
/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/backends/CUT.h>

#include <rp/fts/montecarlo/RealMove.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RealMove<1,CUT>;
      template class RealMove<2,CUT>;
      template class RealMove<3,CUT>;
   }
}
