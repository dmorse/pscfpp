/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryStructureFactor_c.h"

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/complex.h>

#include <rp/fts/analyzer/BinaryStructureFactorBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactor<D,CPT>::BinaryStructureFactor(
                                  Simulator<D,CPT>& simulator,
                                  System<D,CPT>& system)
    : BinaryStructureFactorBase<D,CPT>(simulator, system)
   {}

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D,CPT>::setup()
   {
      allocate();
      WaveList<D,CPT> const & waveList = AnalyzerT::system().waveList();
      findWaveBunches(waveList.kSq(), waveList.implicitInverse());
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D,CPT>::sample(long iStep)
   {
      if (AnalyzerT::isAtInterval(iStep)) {
         computeW();
         computeS(wk_);
      }
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryStructureFactorBase<1,CPT>;
      template class BinaryStructureFactorBase<2,CPT>;
      template class BinaryStructureFactorBase<3,CPT>;
      template class BinaryStructureFactor<1,CPT>;
      template class BinaryStructureFactor<2,CPT>;
      template class BinaryStructureFactor<3,CPT>;
   }
}
