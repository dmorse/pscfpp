/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMove.h"

#include <rpc/fts/montecarlo/McMove.h>
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/montecarlo/ForceBiasMoveBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D, CppTp<D> >::ForceBiasMove(
                              McSimulator<D, CppTp<D> >& simulator)
    : ForceBiasMoveBase<D, CppTp<D> > (simulator)
   {}

   /*
   * Compute force bias field for use in Metropolis acceptance test.
   */
   template<int D>
   void ForceBiasMove<D, CppTp<D> >::computeForceBias(
                               RField<D, CppTp<D> >& result,
                               RField<D, CppTp<D> > const & di,
                               RField<D, CppTp<D> > const & df,
                               RField<D, CppTp<D> > const & dwc,
                               double mobility)
   {
      const int n
         = McMove<D, CppTp<D> >::system().domain().mesh().size();
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
      template class ForceBiasMoveBase<1, CppTp<1> >;
      template class ForceBiasMoveBase<2, CppTp<2> >;
      template class ForceBiasMoveBase<3, CppTp<3> >;
      template class ForceBiasMove<1, CppTp<1> >;
      template class ForceBiasMove<2, CppTp<2> >;
      template class ForceBiasMove<3, CppTp<3> >;
   }
}
