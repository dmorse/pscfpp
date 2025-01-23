/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOpFts.h"
#include <pscf/cuda/ThreadGrid.h>

namespace Pscf {
namespace Rpg {
namespace VecOpFts {

   // CUDA kernels: 
   // (defined in anonymous namespace, used only in this file)

   namespace {

      // Rescale array a from [0,1] to [-b, b]
      __global__ void _mcftsScale(cudaReal* a, cudaReal const b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = a[i] * 2 * b - b;
         }
      }

      // Add array b to real part of a and array c to imaginary part of a
      __global__ void _fourierMove(cudaComplex* a, cudaReal const * b, 
                                   cudaReal const * c, const int n) {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i];
            a[i].y += c[i];
         }
      }

      // Compute d field (functional derivative of H[w])
      __global__ void _computeDField(cudaReal* d, cudaReal const * Wc, 
                                     cudaReal const * Cc, cudaReal const a, 
                                     cudaReal const b, cudaReal const s, 
                                     const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            d[i] = a * (b * (Wc[i] - s) + Cc[i]);
         }
      }

      // Compute force bias
      __global__ void _computeForceBias(cudaReal* result, cudaReal const * di, 
                                        cudaReal const * df, 
                                        cudaReal const * dwc, 
                                        cudaReal mobility, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            result[i] = 0.5 * (di[i] + df[i]) * 
                        (dwc[i] + mobility * (0.5 * (di[i] - df[i])));
         }
      }

   }

   // Kernel wrappers:

   /*
   * Rescale array a from [0,1] to [-b, b], GPU kernel wrapper.
   */
   void mcftsScale(DeviceArray<cudaReal>& a, cudaReal const b)
   {
      const int n = a.capacity();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mcftsScale<<<nBlocks, nThreads>>>(a.cArray(), b, n);
   }

   /*
   * Add array b to real part of a and array c to imaginary part of a
   */
   void fourierMove(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _fourierMove<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), 
                                          c.cArray(), n);
   }

   /*
   * Compute d field (functional derivative of H[w])
   */
   void computeDField(DeviceArray<cudaReal>& d, 
                      DeviceArray<cudaReal> const & Wc, 
                      DeviceArray<cudaReal> const & Cc, 
                      cudaReal const a, cudaReal const b, cudaReal const s)
   {
      const int n = d.capacity();
      UTIL_CHECK(Wc.capacity() == n);
      UTIL_CHECK(Cc.capacity() == n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _computeDField<<<nBlocks, nThreads>>>(d.cArray(), Wc.cArray(), 
                                            Cc.cArray(), a, b, s, n);
   }

   /*
   * Compute force bias
   */
   void computeForceBias(DeviceArray<cudaReal>& result, 
                         DeviceArray<cudaReal> const & di, 
                         DeviceArray<cudaReal> const & df, 
                         DeviceArray<cudaReal> const & dwc, cudaReal mobility)
   {
      const int n = result.capacity();
      UTIL_CHECK(di.capacity() == n);
      UTIL_CHECK(df.capacity() == n);
      UTIL_CHECK(dwc.capacity() == n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _computeForceBias<<<nBlocks, nThreads>>>(result.cArray(), di.cArray(), 
                                               df.cArray(), dwc.cArray(), 
                                               mobility, n);
   }

}
}
}