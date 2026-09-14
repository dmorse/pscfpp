/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourthOrderParameter_u.h"

#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/Reduce.h>
#include <pscf/backend/cuda/HostArray.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cpp/VecOp.h>

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
      // Precondition - check allocation
      int kSize = Base::kSize_;
      UTIL_CHECK(Base::prefactor_.capacity() == kSize);

      // Initialize host array
      HostArray<cudaReal,CUT> prefactor_h;
      prefactor_h.associate(Base::prefactor_);
      UTIL_CHECK(prefactor_h.capacity() == Base::kSize_);

      // Perform computation on host
      Base::computePrefactor(prefactor_h);

      // Copy to device and dissociate host array
      Base::prefactor_ = prefactor_h;
      prefactor_h.dissociate();
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
