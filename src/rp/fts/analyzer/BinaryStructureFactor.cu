/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>
#include <pscf/backend/cuda/ConstHostArray.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/complex.h>

#include <rp/fts/analyzer/BinaryStructureFactor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryStructureFactor<1,CUT>;
      template class BinaryStructureFactor<2,CUT>;
      template class BinaryStructureFactor<3,CUT>;
   }
}
