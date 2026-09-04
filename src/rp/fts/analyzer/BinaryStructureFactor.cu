/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryStructureFactor_u.h"

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

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D,CUT>::setup()
   {
      allocate();
      UTIL_CHECK(wk_.isAllocated());
      if (!wkHost_.isAllocated()) {
         wkHost_.allocate(wk_.capacity());
      }

      WaveList<D,CUT> const & waveList = AnalyzerT::system().waveList();
      HostDArray<double> kSq = waveList.kSq();
      HostDArray<bool> implicit = waveList.implicitInverse();
      Base::findWaveBunches(kSq, implicit);
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D,CUT>::sample(long iStep)
   {
      if (AnalyzerT::isAtInterval(iStep)) {
         Base::computeW();
         wkHost_ = wk_;  // Copy wk_ from device to host
         Base::computeS(wkHost_);
      }
   }

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
