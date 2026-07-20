/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryStructureFactor.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>

#include <prdc/field/cuda/FFT.h>
#include <prdc/field/cuda/WaveList.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/complex.h>

#include <rp/fts/analyzer/BinaryStructureFactorBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactor<D, CudaTp<D> >::BinaryStructureFactor(
              Simulator<D, CudaTp<D> >& simulator,
              System<D, CudaTp<D> >& system)
    : BinaryStructureFactorBase< D, CudaTp<D> >(simulator, system)
   {}

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D, CudaTp<D> >::setup()
   {
      allocate();
      UTIL_CHECK(wk_.isAllocated());
      if (!wkHost_.isAllocated()) {
         wkHost_.allocate(wk_.capacity());
      }

      Cuda::WaveList<D> const & waveList = AnalyzerT::system().waveList();
      HostDArray<double> kSq = waveList.kSq();
      HostDArray<bool> implicit = waveList.implicitInverse();
      Base::findWaveBunches(kSq, implicit);
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D, CudaTp<D> >::sample(long iStep)
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
      template class BinaryStructureFactorBase<1, CudaTp<1> >;
      template class BinaryStructureFactorBase<2, CudaTp<2> >;
      template class BinaryStructureFactorBase<3, CudaTp<3> >;
      template class BinaryStructureFactor<1, CudaTp<1> >;
      template class BinaryStructureFactor<2, CudaTp<2> >;
      template class BinaryStructureFactor<3, CudaTp<3> >;
   }
}
