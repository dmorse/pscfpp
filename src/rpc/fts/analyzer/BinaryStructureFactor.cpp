/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryStructureFactor.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/WaveList.h>

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
   BinaryStructureFactor<D, Rpc::Types<D> >::BinaryStructureFactor(
                                  Simulator<D, Rpc::Types<D> >& simulator,
                                  System<D, Rpc::Types<D> >& system)
    : BinaryStructureFactorBase< D, Rpc::Types<D> >(simulator, system)
   {}

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D, Rpc::Types<D> >::setup()
   {
      allocate();
      Cpu::WaveList<D> const & waveList = AnalyzerT::system().waveList();
      findWaveBunches(waveList.kSq(), waveList.implicitInverse());
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D, Rpc::Types<D> >::sample(long iStep)
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
      template class BinaryStructureFactorBase<1, Rpc::Types<1> >;
      template class BinaryStructureFactorBase<2, Rpc::Types<2> >;
      template class BinaryStructureFactorBase<3, Rpc::Types<3> >;
      template class BinaryStructureFactor<1, Rpc::Types<1> >;
      template class BinaryStructureFactor<2, Rpc::Types<2> >;
      template class BinaryStructureFactor<3, Rpc::Types<3> >;
   }
}
