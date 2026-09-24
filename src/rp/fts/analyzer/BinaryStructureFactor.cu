/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryStructureFactor_u.h"

#include <pscf/backend/cuda/ConstHostArray.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/complex.h>

#include <rp/fts/analyzer/BinaryStructureFactorBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactor<D,CUT>::BinaryStructureFactor(
              Simulator<D,CUT>& simulator,
              System<D,CUT>& system)
    : BinaryStructureFactorBase<D,CUT>(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryStructureFactorBase<1,CUT>;
      template class BinaryStructureFactorBase<2,CUT>;
      template class BinaryStructureFactorBase<3,CUT>;
      template class BinaryStructureFactor<1,CUT>;
      template class BinaryStructureFactor<2,CUT>;
      template class BinaryStructureFactor<3,CUT>;
   }
}
