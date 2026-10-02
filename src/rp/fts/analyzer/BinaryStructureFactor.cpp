/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <pscf/backend/cpp/ConstHostArray.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/VecOpCx.h>
#include <pscf/backend/cpp/complex.h>

#include <rp/fts/analyzer/BinaryStructureFactor.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryStructureFactor<1,CPT>;
      template class BinaryStructureFactor<2,CPT>;
      template class BinaryStructureFactor<3,CPT>;
   }
}
