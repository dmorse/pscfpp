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

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/WaveList.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/complex.h>

#include <rp/fts/analyzer/BinaryStructureFactor.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   BinaryStructureFactor<D>::BinaryStructureFactor(
                                  Simulator<D>& simulator,
                                  Rp::System<D, Rpg::Types<D> >& system)
    : Rp::BinaryStructureFactor< D, Types<D> >(simulator, system)
   {}

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D>::setup()
   {
      allocate();
      UTIL_CHECK(wk_.isAllocated());
      if (!wkHost_.isAllocated()) {
         wkHost_.allocate(wk_.capacity());
      }

      Cuda::WaveList<D> const & waveList = AnalyzerT::system().waveList();
      HostDArray<double> kSq = waveList.kSq();
      HostDArray<bool> implicit = waveList.implicitInverse();
      findWaveBunches(kSq, implicit);
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D>::sample(long iStep)
   {
      if (AnalyzerT::isAtInterval(iStep)) {
         computeW();
         wkHost_ = wk_;  // Copy wk_ from device to host
         computeS(wkHost_);
      }
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryStructureFactor<1, Rpg::Types<1> >;
      template class BinaryStructureFactor<2, Rpg::Types<2> >;
      template class BinaryStructureFactor<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class BinaryStructureFactor<1>;
      template class BinaryStructureFactor<2>;
      template class BinaryStructureFactor<3>;
   }
}
