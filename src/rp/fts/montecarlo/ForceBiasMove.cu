/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaTypes.h>

#include <rp/fts/montecarlo/ForceBiasMove_u.h>
#include <rp/fts/montecarlo/ForceBiasMoveBase.tpp>

// Anonymous namespace for CUDA kernel
namespace Pscf {
namespace {

   // Compute force bias
   __global__
   void _computeForceBias(cudaReal* result,
                          cudaReal const * di,
                          cudaReal const * df,
                          cudaReal const * dwc,
                          double mobility,
                          const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         result[i] = 0.5 * (di[i] + df[i]) *
                     (dwc[i] + mobility * (0.5 * (di[i] - df[i])));
      }
   }

}
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D,CUT>::ForceBiasMove(
                             McSimulator<D,CUT>& simulator)
    : ForceBiasMoveBase<D,CUT> (simulator)
   {}

   /*
   * Compute force bias field for use in Metropolis acceptance test.
   */
   template<int D>
   void ForceBiasMove<D,CUT>::computeForceBias(
                               RField<D,CUT>& result,
                               RField<D,CUT> const & di,
                               RField<D,CUT> const & df,
                               RField<D,CUT> const & dwc,
                               double mobility)
   {
      const int n = result.capacity();
      UTIL_CHECK(di.capacity() == n);
      UTIL_CHECK(df.capacity() == n);
      UTIL_CHECK(dwc.capacity() == n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _computeForceBias<<<nBlocks, nThreads>>>(
                        result.cArray(), di.cArray(),
                        df.cArray(), dwc.cArray(),
                        mobility, n);
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ForceBiasMoveBase<1,CUT>;
      template class ForceBiasMoveBase<2,CUT>;
      template class ForceBiasMoveBase<3,CUT>;
      template class ForceBiasMove<1,CUT>;
      template class ForceBiasMove<2,CUT>;
      template class ForceBiasMove<3,CUT>;
   }
}
