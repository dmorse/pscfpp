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
   BinaryStructureFactor<D, CppTp<D> >::BinaryStructureFactor(
                                  Simulator<D, CppTp<D> >& simulator,
                                  System<D, CppTp<D> >& system)
    : BinaryStructureFactorBase< D, CppTp<D> >(simulator, system)
   {}

   /*
   * Setup before entering main loop.
   */
   template <int D>
   void BinaryStructureFactor<D, CppTp<D> >::setup()
   {
      allocate();
      WaveList<D, CppTp<D> > const & waveList = AnalyzerT::system().waveList();
      findWaveBunches(waveList.kSq(), waveList.implicitInverse());
   }

   /*
   * Compute structure factors for all wavevectors and bunches.
   */
   template <int D>
   void BinaryStructureFactor<D, CppTp<D> >::sample(long iStep)
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
      template class BinaryStructureFactorBase<1, CppTp<1> >;
      template class BinaryStructureFactorBase<2, CppTp<2> >;
      template class BinaryStructureFactorBase<3, CppTp<3> >;
      template class BinaryStructureFactor<1, CppTp<1> >;
      template class BinaryStructureFactor<2, CppTp<2> >;
      template class BinaryStructureFactor<3, CppTp<3> >;
   }
}
