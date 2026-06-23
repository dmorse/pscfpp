/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourthOrderParameter.h"

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <prdc/cuda/FFT.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaTypes.h>
#include <pscf/cpu/VecOp.h>

#include <rp/fts/analyzer/FourthOrderParameter.tpp>

namespace Pscf {
namespace Rpg {

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D>::FourthOrderParameter(
                                   Simulator<D>& simulator,
                                   Rp::System<D, Rpg::Types<D> >& system)
    : Base(simulator, system)
   {}

   /*
   * Initialize Base::prefactor_ protected member variable.
   */
   template <int D>
   void FourthOrderParameter<D>::computePrefactor()
   {
      // Allocate CPU host array
      HostDArray<cudaReal> prefactor_h(Base::kSize_);
      VecOp::eqS(prefactor_h, 0.0);

      // Perform computation on host
      Base::computePrefactor(prefactor_h);

      // Copy from from cpu(host) to gpu(device)
      Base::prefactor_ = prefactor_h;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FourthOrderParameter<1, Rpg::Types<1> >;
      template class FourthOrderParameter<2, Rpg::Types<2> >;
      template class FourthOrderParameter<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class FourthOrderParameter<1>;
      template class FourthOrderParameter<2>;
      template class FourthOrderParameter<3>;
   }
}
