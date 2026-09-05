/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>
#include <pscf/backend/cpp/CpuVecRandom.h>

#include "ForceBiasMove_c.h"
#include <rp/fts/montecarlo/ForceBiasMoveBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D,CPT>::ForceBiasMove(
                              McSimulator<D,CPT>& simulator)
    : ForceBiasMoveBase<D,CPT> (simulator)
   {}

   /*
   * Compute force bias field for use in Metropolis acceptance test.
   */
   template<int D>
   void ForceBiasMove<D,CPT>::computeForceBias(
                               RField<D,CPT>& result,
                               RField<D,CPT> const & di,
                               RField<D,CPT> const & df,
                               RField<D,CPT> const & dwc,
                               double mobility)
   {
      const int n
         = McMove<D,CPT>::system().domain().mesh().size();
      UTIL_CHECK(result.capacity() == n);
      UTIL_CHECK(di.capacity() == n);
      UTIL_CHECK(df.capacity() == n);
      UTIL_CHECK(dwc.capacity() == n);

      double dp, dm;
      for (int k = 0; k < n; ++k) {
         dp = 0.5*(di[k] + df[k]);
         dm = 0.5*(di[k] - df[k]);
         result[k] = dp*( dwc[k] + mobility*dm );
      }
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ForceBiasMoveBase<1,CPT>;
      template class ForceBiasMoveBase<2,CPT>;
      template class ForceBiasMoveBase<3,CPT>;
      template class ForceBiasMove<1,CPT>;
      template class ForceBiasMove<2,CPT>;
      template class ForceBiasMove<3,CPT>;
   }
}
