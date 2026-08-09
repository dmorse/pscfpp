/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/backends/CPT.h>

#include <rp/fts/compressor/LrAmCompressor.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LrAmCompressor<1,CPT>;
      template class LrAmCompressor<2,CPT>;
      template class LrAmCompressor<3,CPT>;
   }
}
