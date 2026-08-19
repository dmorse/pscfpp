/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourthOrderParameter_u.h"

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaTypes.h>
#include <pscf/cpu/VecOp.h>

#include <rp/fts/analyzer/FourthOrderParameterBase.tpp>

namespace Pscf {
namespace Rp {

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D,CUT>::FourthOrderParameter(
                                   Rp::Simulator<D,CUT>& simulator,
                                   Rp::System<D,CUT>& system)
    : Base(simulator, system)
   {}

   /*
   * Initialize Base::prefactor_ protected member variable.
   */
   template <int D>
   void FourthOrderParameter<D,CUT>::computePrefactor()
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
      template class FourthOrderParameterBase<1,CUT>;
      template class FourthOrderParameterBase<2,CUT>;
      template class FourthOrderParameterBase<3,CUT>;
      template class FourthOrderParameter<1,CUT>;
      template class FourthOrderParameter<2,CUT>;
      template class FourthOrderParameter<3,CUT>;
   }
}
