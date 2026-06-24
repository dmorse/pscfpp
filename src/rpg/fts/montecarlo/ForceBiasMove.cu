/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMove.h"
#include "McMove.h"
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>

#include <rp/fts/montecarlo/ForceBiasMove.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   ForceBiasMove<D>::ForceBiasMove(Rp::McSimulator<D, Rpg::Types<D> >& simulator)
    : Rp::ForceBiasMove<D, Types<D> > (simulator)
   {}

   // Anonymous namespace for CUDA kernel
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

   /*
   * Compute force bias field for use in Metropolis acceptance test.
   */
   template<int D>
   void ForceBiasMove<D>::computeForceBias(
                               RField<D>& result,
                               RField<D> const & di,
                               RField<D> const & df,
                               RField<D> const & dwc,
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
      template class ForceBiasMove<1, Rpg::Types<1> >;
      template class ForceBiasMove<2, Rpg::Types<2> >;
      template class ForceBiasMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class ForceBiasMove<1>;
      template class ForceBiasMove<2>;
      template class ForceBiasMove<3>;
   }
}
