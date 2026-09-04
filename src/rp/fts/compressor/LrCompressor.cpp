/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOpCx.h>
#include <pscf/backend/cpp/Reduce.h>
#include <pscf/backend/cpp/CPT.h>

#include <rp/fts/compressor/LrCompressor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrCompressor<1,CPT>;
      template class LrCompressor<2,CPT>;
      template class LrCompressor<3,CPT>;
   }
}
