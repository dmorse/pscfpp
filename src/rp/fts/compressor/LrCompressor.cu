/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>
#include <pscf/backend/cuda/CUT.h>

#include <rp/fts/compressor/LrCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrCompressor<1,CUT>;
      template class LrCompressor<2,CUT>;
      template class LrCompressor<3,CUT>;
   }
}
