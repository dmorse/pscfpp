/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Reduce.h"
#include "ThreadGrid.h"
#include "HostDArray.h"
#include <cmath>

namespace Pscf {
namespace Reduce {

/*
* Compute sum of array elements (GPU kernel wrapper).
*/
__host__ cudaReal sum(DeviceArray<cudaReal> const & in) 
{
   UTIL_CHECK(in.isAllocated());
   int n = in.capacity();

   // Set up temporary device arrays for storing reduced data
   DeviceArray<cudaReal> temp1, temp2;
   
   int i = 0;

   // Perform parallel reduction on GPU repeatedly until n < 1e5
   while (n >= 1e5) {
      // Establish GPU resources for this parallel reduction. Divided by 
      // two because of the global memory load in the kernel performing 
      // the first level of reduction!
      int nBlocks, nThreads;
      int halvedSize = ceil((float)n/2);

      ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);
      // Note: setThreadsLogical ensures that nThreads is a power of 2

      if (nThreads < 64) {
         // Thread blocks too small. Manually set nThreads to 64
         ThreadGrid::setThreadsPerBlock(64);
         ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for parallel reduction.\n"
                     << "Setting nThreads equal to 64." << std::endl;
      }

      // Warp size must be 32
      UTIL_CHECK(ThreadGrid::warpSize() == 32);

      // Perform parallel reduction
      if (i == 0) { // first reduction, use input array
         temp1.allocate(nBlocks);
         _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), in.cArray(), n);
      } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
         temp2.allocate(nBlocks);
         _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp2.cArray(), temp1.cArray(), n);
         temp1.deallocate();
      } else {                 // i is even: reduce temp2, store in temp1
         temp1.allocate(nBlocks);
         _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), temp2.cArray(), n);
         temp2.deallocate();
      }

      n = nBlocks;
      i += 1;
   }

   // Transfer the partially reduced sum to the host
   HostDArray<cudaReal> temp_h;
   if (i == 0) {
      temp_h = in;
   } else if (i % 2 == 1) {
      temp_h = temp1;
   } else {
      temp_h = temp2;
   }

   if (n == 1) {
      return temp_h[0];
   } else {
      // Sum up elements of temp_h to get final result.
      // Use Kahan summation to reduce accumulation of error
      cudaReal sum = 0.0, tempVal, tempSum;
      cudaReal err = 0.0;
      for (int i = 0; i < n; ++i) {
         tempVal = temp_h[i] - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;  
      }
      return sum;
   }
}

/*
* Get maximum of array elements (GPU kernel wrapper).
*/
__host__ cudaReal max(DeviceArray<cudaReal> const & in)
{
   UTIL_CHECK(in.isAllocated());
   int n = in.capacity();

   // Set up temporary device arrays for storing reduced data
   DeviceArray<cudaReal> temp1, temp2;
   
   int i = 0;

   // Perform parallel reduction on GPU repeatedly until n < 1e5
   while (n >= 1e5) {
      // Establish GPU resources for this parallel reduction. Divided by 
      // two because of the global memory load in the kernel performing 
      // the first level of reduction!
      int nBlocks, nThreads;
      int halvedSize = ceil((float)n/2);

      ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);
      // Note: setThreadsLogical ensures that nThreads is a power of 2

      if (nThreads < 64) {
         // Thread blocks too small. Manually set nThreads to 64
         ThreadGrid::setThreadsPerBlock(64);
         ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for parallel reduction.\n"
                     << "Setting nThreads equal to 64." << std::endl;
      }

      // Warp size must be 32
      UTIL_CHECK(ThreadGrid::warpSize() == 32);

      // Perform parallel reduction
      if (i == 0) { // first reduction, use input array
         temp1.allocate(nBlocks);
         _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), in.cArray(), n);
      } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
         temp2.allocate(nBlocks);
         _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp2.cArray(), temp1.cArray(), n);
         temp1.deallocate();
      } else {                 // i is even: reduce temp2, store in temp1
         temp1.allocate(nBlocks);
         _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), temp2.cArray(), n);
         temp2.deallocate();
      }

      n = nBlocks;
      i += 1;
   }

   // Transfer the partially reduced sum to the host
   HostDArray<cudaReal> temp_h;
   if (i == 0) {
      temp_h = in;
   } else if (i % 2 == 1) {
      temp_h = temp1;
   } else {
      temp_h = temp2;
   }

   cudaReal max = temp_h[0];
   for (int i = 1; i < n; i++) {
      if (temp_h[i] > max) max = temp_h[i];
   }
   return max;
}

/*
* Get maximum absolute magnitude of array elements (GPU kernel wrapper).
*/
__host__ cudaReal maxAbs(DeviceArray<cudaReal> const & in)
{
   UTIL_CHECK(in.isAllocated());
   int n = in.capacity();

   // Set up temporary device arrays for storing reduced data
   DeviceArray<cudaReal> temp1, temp2;
   
   int i = 0;

   // Perform parallel reduction on GPU repeatedly until n < 1e5
   //   (note: for this wrapper, we always call the kernel at least once,
   //    even if n < 1e5, so that the part done on the CPU is always just
   //    comparing array element size, without needing fabs().)
   while (n >= 1e5 || i == 0) {
      // Establish GPU resources for this parallel reduction. Divided by 
      // two because of the global memory load in the kernel performing 
      // the first level of reduction!
      int nBlocks, nThreads;
      int halvedSize = ceil((float)n/2);

      ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);
      // Note: setThreadsLogical ensures that nThreads is a power of 2

      if (nThreads < 64) {
         // Thread blocks too small. Manually set nThreads to 64
         ThreadGrid::setThreadsPerBlock(64);
         ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for parallel reduction.\n"
                     << "Setting nThreads equal to 64." << std::endl;
      }

      // Warp size must be 32
      UTIL_CHECK(ThreadGrid::warpSize() == 32);

      // Perform parallel reduction
      if (i == 0) { // first reduction, use input array
         temp1.allocate(nBlocks);
         _maxAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), in.cArray(), n);
      } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
         temp2.allocate(nBlocks);
         _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp2.cArray(), temp1.cArray(), n);
         temp1.deallocate();
      } else {                 // i is even: reduce temp2, store in temp1
         temp1.allocate(nBlocks);
         _max<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), temp2.cArray(), n);
         temp2.deallocate();
      }

      n = nBlocks;
      i += 1;
   }

   // Transfer the partially reduced sum to the host
   HostDArray<cudaReal> temp_h;
   if (i == 0) {
      temp_h = in;
   } else if (i % 2 == 1) {
      temp_h = temp1;
   } else {
      temp_h = temp2;
   }

   cudaReal max = temp_h[0];
   for (int i = 1; i < n; i++) {
      if (temp_h[i] > max) max = temp_h[i];
   }
   return max;
}

/*
* Get minimum of array elements (GPU kernel wrapper).
*/
__host__ cudaReal min(DeviceArray<cudaReal> const & in)
{
   UTIL_CHECK(in.isAllocated());
   int n = in.capacity();

   // Set up temporary device arrays for storing reduced data
   DeviceArray<cudaReal> temp1, temp2;
   
   int i = 0;

   // Perform parallel reduction on GPU repeatedly until n < 1e5
   while (n >= 1e5) {
      // Establish GPU resources for this parallel reduction. Divided by 
      // two because of the global memory load in the kernel performing 
      // the first level of reduction!
      int nBlocks, nThreads;
      int halvedSize = ceil((float)n/2);

      ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);
      // Note: setThreadsLogical ensures that nThreads is a power of 2

      if (nThreads < 64) {
         // Thread blocks too small. Manually set nThreads to 64
         ThreadGrid::setThreadsPerBlock(64);
         ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for parallel reduction.\n"
                     << "Setting nThreads equal to 64." << std::endl;
      }

      // Warp size must be 32
      UTIL_CHECK(ThreadGrid::warpSize() == 32);

      // Perform parallel reduction
      if (i == 0) { // first reduction, use input array
         temp1.allocate(nBlocks);
         _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), in.cArray(), n);
      } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
         temp2.allocate(nBlocks);
         _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp2.cArray(), temp1.cArray(), n);
         temp1.deallocate();
      } else {                 // i is even: reduce temp2, store in temp1
         temp1.allocate(nBlocks);
         _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), temp2.cArray(), n);
         temp2.deallocate();
      }

      n = nBlocks;
      i += 1;
   }

   // Transfer the partially reduced sum to the host
   HostDArray<cudaReal> temp_h;
   if (i == 0) {
      temp_h = in;
   } else if (i % 2 == 1) {
      temp_h = temp1;
   } else {
      temp_h = temp2;
   }

   cudaReal min = temp_h[0];
   for (int i = 1; i < n; i++) {
      if (temp_h[i] < min) min = temp_h[i];
   }
   return min;
}

/*
* Get minimum absolute magnitude of array elements (GPU kernel wrapper).
*/
__host__ cudaReal minAbs(DeviceArray<cudaReal> const & in)
{
   UTIL_CHECK(in.isAllocated());
   int n = in.capacity();

   // Set up temporary device arrays for storing reduced data
   DeviceArray<cudaReal> temp1, temp2;
   
   int i = 0;

   // Perform parallel reduction on GPU repeatedly until n < 1e5
   //   (note: for this wrapper, we always call the kernel at least once,
   //    even if n < 1e5, so that the part done on the CPU is always just
   //    comparing array element size, without needing fabs().)
   while (n >= 1e5 || i == 0) {
      // Establish GPU resources for this parallel reduction. Divided by 
      // two because of the global memory load in the kernel performing 
      // the first level of reduction!
      int nBlocks, nThreads;
      int halvedSize = ceil((float)n/2);

      ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);
      // Note: setThreadsLogical ensures that nThreads is a power of 2

      if (nThreads < 64) {
         // Thread blocks too small. Manually set nThreads to 64
         ThreadGrid::setThreadsPerBlock(64);
         ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for parallel reduction.\n"
                     << "Setting nThreads equal to 64." << std::endl;
      }

      // Warp size must be 32
      UTIL_CHECK(ThreadGrid::warpSize() == 32);

      // Perform parallel reduction
      if (i == 0) { // first reduction, use input array
         temp1.allocate(nBlocks);
         _minAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), in.cArray(), n);
      } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
         temp2.allocate(nBlocks);
         _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp2.cArray(), temp1.cArray(), n);
         temp1.deallocate();
      } else {                 // i is even: reduce temp2, store in temp1
         temp1.allocate(nBlocks);
         _min<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), temp2.cArray(), n);
         temp2.deallocate();
      }

      n = nBlocks;
      i += 1;
   }

   // Transfer the partially reduced sum to the host
   HostDArray<cudaReal> temp_h;
   if (i == 0) {
      temp_h = in;
   } else if (i % 2 == 1) {
      temp_h = temp1;
   } else {
      temp_h = temp2;
   }

   cudaReal min = temp_h[0];
   for (int i = 1; i < n; i++) {
      if (temp_h[i] < min) min = temp_h[i];
   }
   return min;
}

/*
* Compute inner product of two real arrays (GPU kernel wrapper).
*/
__host__ cudaReal innerProduct(DeviceArray<cudaReal> const & a,  
                               DeviceArray<cudaReal> const & b)
{
   UTIL_CHECK(a.isAllocated());
   UTIL_CHECK(b.isAllocated());
   UTIL_CHECK(a.capacity() == b.capacity());
   int n = a.capacity();

   // Set up temporary device arrays for storing reduced data
   DeviceArray<cudaReal> temp1, temp2;
   
   int i = 0;

   // Perform parallel reduction on GPU repeatedly until n < 1e5
   //   (note: for this wrapper, we always call the kernel at least once,
   //    even if n < 1e5, so that the part done on the CPU is always just
   //    adding up array elements.)
   while (n >= 1e5 || i == 0) {
      // Establish GPU resources for this parallel reduction. Divided by 
      // two because of the global memory load in the kernel performing 
      // the first level of reduction!
      int nBlocks, nThreads;
      int halvedSize = ceil((float)n/2);

      ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);
      // Note: setThreadsLogical ensures that nThreads is a power of 2

      if (nThreads < 64) {
         // Thread blocks too small. Manually set nThreads to 64
         ThreadGrid::setThreadsPerBlock(64);
         ThreadGrid::setThreadsLogical(halvedSize,nBlocks,nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for parallel reduction.\n"
                     << "Setting nThreads equal to 64." << std::endl;
      }

      // Warp size must be 32
      UTIL_CHECK(ThreadGrid::warpSize() == 32);

      // Perform parallel reduction

      // (note: only the first kernel call uses _innerProduct. After
      //  that, we use _sum to reduce the array output by _innerProduct.)

      if (i == 0) { // first reduction, use input arrays
         temp1.allocate(nBlocks);
         _innerProduct<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                              (temp1.cArray(), a.cArray(), b.cArray(), n);
      } else if (i % 2 == 1) { // i is odd: reduce temp1, store in temp2
         temp2.allocate(nBlocks);
         _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp2.cArray(), temp1.cArray(), n);
         temp1.deallocate();
      } else {                 // i is even: reduce temp2, store in temp1
         temp1.allocate(nBlocks);
         _sum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (temp1.cArray(), temp2.cArray(), n);
         temp2.deallocate();
      }

      n = nBlocks;
      i += 1;
   }

   // Transfer the partially reduced sum to the host
   HostDArray<cudaReal> temp_h;
   if (i % 2 == 1) {
      temp_h = temp1;
   } else {
      temp_h = temp2;
   }

   if (n == 1) {
      return temp_h[0];
   } else {
      // Sum up elements of temp_h to get final result.
      // Use Kahan summation to reduce accumulation of error
      cudaReal sum = 0.0, tempVal, tempSum;
      cudaReal err = 0.0;
      for (int i = 0; i < n; ++i) {
         tempVal = temp_h[i] - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;  
      }
      return sum;
   }
}

/*
* Compute sum of array elements (GPU kernel).
*/
__global__ void _sum(cudaReal* sum, const cudaReal* in, int n)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tId = threadIdx.x;
   int bId = blockIdx.x;
   int bDim = blockDim.x;
   int idx = bId * (bDim*2) + tId;

   // Shared memory holding area
   extern __shared__ cudaReal sData[];

   // Global memory load and first operation
   if (idx < n) {
      sData[tId] = in[idx];
      if (idx + bDim < n) {
         sData[tId] += in[idx+bDim];
      }
   } else {
      // idx > n. Set value to 0.0 to fully populate sData without 
      // contributing to the sum
      sData[tId] = 0.0;
   }
   
   // wait for all threads to finish
   __syncthreads();

   // Make reductions across the block of data, each thread handling 
   // one reduction across two data points with strided indices before 
   // syncing with each other and then making further reductions.
   for (int stride = bDim / 2; stride > 32; stride /= 2) {
      if (tId < stride) {
         sData[tId] += sData[tId+stride];
      }
      __syncthreads();
   }

   // Unwrap last warp (stride == 32)
   if (tId < 32) {
      _warpSum(sData, tId); // defined at bottom of this file
   }

   // Store the output of the threads in this block
   if (tId == 0) {
      sum[bId] = sData[0];
   }
}

/*
* Get maximum of array elements (GPU kernel).
*/
__global__ void _max(cudaReal* max, const cudaReal* in, int n)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tId = threadIdx.x;
   int bId = blockIdx.x;
   int bDim = blockDim.x;
   int idx = bId * (bDim*2) + tId;
   
   // Shared memory holding area
   extern __shared__ cudaReal sData[];

   // Global memory load and first operation
   if (idx < n) {
      sData[tId] = in[idx];
      if (idx + bDim < n) {
         cudaReal in1 = in[idx+bDim];
         sData[tId] = (sData[tId] > in1) ? sData[tId] : in1;
      }
   } else {
      // idx > n. Set value to in[idx-n], an earlier value in the 
      // array, to fully populate sData without altering the result
      sData[tId] = in[idx-n];
   }
   
   // wait for all threads to finish
   __syncthreads();

   // Make reductions across the block of data, each thread handling 
   // one reduction across two data points with strided indices before 
   // syncing with each other and then making further reductions.
   for (int stride = bDim/2; stride > 32; stride/=2) {
      if (tId < stride) {
         if (sData[tId+stride] > sData[tId]) {
            sData[tId] = sData[tId+stride];
         }
      }
      __syncthreads();
   }

   // Unwrap last warp (stride == 32)
   if (tId < 32) {
      _warpMax(sData, tId); // defined at bottom of this file
   }

   // Store the output of the threads in this block
   if (tId == 0) {
      max[bId] = sData[0];
   }
}

/*
* Get maximum absolute magnitude of array elements (GPU kernel).
*/
__global__ void _maxAbs(cudaReal* max, const cudaReal* in, int n) 
{
   // number of blocks cut in two to avoid inactive initial threads
   int tId = threadIdx.x;
   int bId = blockIdx.x;
   int bDim = blockDim.x;
   int idx = bId * (bDim*2) + tId;
   
   // Shared memory holding area
   extern __shared__ cudaReal sData[];

   // Global memory load and first operation
   if (idx < n) {
      sData[tId] = fabs(in[idx]);
      if (idx + bDim < n) {
         cudaReal in1 = fabs(in[idx+bDim]);
         sData[tId] = (sData[tId] > in1) ? sData[tId] : in1;
      }
   } else {
      // idx > n. Set value to 0.0 to fully populate sData without 
      // altering the result
      sData[tId] = 0.0;
   }
   
   // wait for all threads to finish
   __syncthreads();

   // Make reductions across the block of data, each thread handling 
   // one reduction across two data points with strided indices before 
   // syncing with each other and then making further reductions.
   for (int stride = bDim/2; stride > 32; stride/=2) {
      if (tId < stride) {
         if (sData[tId+stride] > sData[tId]) {
            sData[tId] = sData[tId+stride];
         }
      }
      __syncthreads();
   }

   // Unwrap last warp (stride == 32)
   if (tId < 32) {
      _warpMax(sData, tId); // defined at bottom of this file
   }

   // Store the output of the threads in this block
   if (tId == 0) {
      max[bId] = sData[0];
   }
}

/*
* Get minimum of array elements (GPU kernel).
*/
__global__ void _min(cudaReal* min, const cudaReal* in, int n)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tId = threadIdx.x;
   int bId = blockIdx.x;
   int bDim = blockDim.x;
   int idx = bId * (bDim*2) + tId;
   
   // Shared memory holding area
   extern __shared__ cudaReal sData[];

   // Global memory load and first operation
   if (idx < n) {
      sData[tId] = in[idx];
      if (idx + bDim < n) {
         cudaReal in1 = in[idx+bDim];
         sData[tId] = (sData[tId] < in1) ? sData[tId] : in1;
      }
   } else {
      // idx > n. Set value to in[idx-n], an earlier value in the 
      // array, to fully populate sData without altering the result
      sData[tId] = in[idx-n];
   }
   
   // wait for all threads to finish
   __syncthreads();

   // Make reductions across the block of data, each thread handling 
   // one reduction across two data points with strided indices before 
   // syncing with each other and then making further reductions.
   for (int stride = bDim/2; stride > 32; stride/=2) {
      if (tId < stride) {
         if (sData[tId+stride] < sData[tId]) {
            sData[tId] = sData[tId+stride];
         }
      }
      __syncthreads();
   }

   // Unwrap last warp (stride == 32)
   if (tId < 32) {
      _warpMin(sData, tId); // defined at bottom of this file
   }

   // Store the output of the threads in this block
   if (tId == 0) {
      min[bId] = sData[0];
   }
}

/*
* Get minimum absolute magnitude of array elements (GPU kernel).
*/
__global__ void _minAbs(cudaReal* min, const cudaReal* in, int n)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tId = threadIdx.x;
   int bId = blockIdx.x;
   int bDim = blockDim.x;
   int idx = bId * (bDim*2) + tId;
   
   // Shared memory holding area
   extern __shared__ cudaReal sData[];

   // Global memory load and first operation
   if (idx < n) {
      sData[tId] = fabs(in[idx]);
      if (idx + bDim < n) {
         cudaReal in1 = fabs(in[idx+bDim]);
         sData[tId] = (sData[tId] < in1) ? sData[tId] : in1;
      }
   } else {
      // idx > n. Set value to fabs(in[idx-n]), an earlier value in the 
      // array, to fully populate sData without altering the result
      sData[tId] = fabs(in[idx-n]);
   }
   
   // wait for all threads to finish
   __syncthreads();

   // Make reductions across the block of data, each thread handling 
   // one reduction across two data points with strided indices before 
   // syncing with each other and then making further reductions.
   for (int stride = bDim/2; stride > 32; stride/=2) {
      if (tId < stride) {
         if (sData[tId+stride] < sData[tId]) {
            sData[tId] = sData[tId+stride];
         }
      }
      __syncthreads();
   }

   // Unwrap last warp (stride == 32)
   if (tId < 32) {
      _warpMin(sData, tId); // defined at bottom of this file
   }

   // Store the output of the threads in this block
   if (tId == 0) {
      min[bId] = sData[0];
   }
}

/*
* Compute inner product of two real arrays (GPU kernel).
*/
__global__ void _innerProduct(cudaReal* ip, const cudaReal* a, 
                              const cudaReal* b, int n)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tId = threadIdx.x;
   int bId = blockIdx.x;
   int bDim = blockDim.x;
   int idx = bId * (bDim*2) + tId;

   // Shared memory holding area
   extern __shared__ cudaReal sData[];

   // Global memory load and first operation
   if (idx < n) {
      sData[tId] = a[idx] * b[idx];
      if (idx + bDim < n) {
         sData[tId] += (a[idx+bDim] * b[idx+bDim]);
      }
   } else { 
      // idx > n. Set value to 0.0 to fully populate sData without 
      // contributing to the sum
      sData[tId] = 0.0;
   }   

   // wait for all threads to finish
   __syncthreads();

   // Make reductions across the block of data, each thread handling 
   // one reduction across two data points with strided indices before 
   // syncing with each other and then making further reductions.
   for (int stride = bDim / 2; stride > 32; stride /= 2) {
      if (tId < stride) {
         sData[tId] += sData[tId+stride];
      }
      __syncthreads();
   }

   // Unwrap last warp (stride == 32)
   if (tId < 32) {
      _warpSum(sData, tId); // defined at bottom of this file
   }

   // Store the output of the threads in this block
   if (tId == 0) {
      ip[bId] = sData[0];
   }
}

/*
* Utility to perform parallel reduction summation within a single warp.
*/
__device__ void _warpSum(volatile cudaReal* sData, int tId)
{
   sData[tId] += sData[tId + 32];
   sData[tId] += sData[tId + 16];
   sData[tId] += sData[tId + 8];
   sData[tId] += sData[tId + 4];
   sData[tId] += sData[tId + 2];
   sData[tId] += sData[tId + 1];
}

/*
* Utility to perform parallel reduction maximization within a single warp.
*/
__device__ void _warpMax(volatile cudaReal* sData, int tId)
{
   if (sData[tId + 32] > sData[tId]) sData[tId] = sData[tId + 32];
   if (sData[tId + 16] > sData[tId]) sData[tId] = sData[tId + 16];
   if (sData[tId + 8] > sData[tId]) sData[tId] = sData[tId + 8];
   if (sData[tId + 4] > sData[tId]) sData[tId] = sData[tId + 4];
   if (sData[tId + 2] > sData[tId]) sData[tId] = sData[tId + 2];
   if (sData[tId + 1] > sData[tId]) sData[tId] = sData[tId + 1];
}

/*
* Utility to perform parallel reduction minimization within a single warp.
*/
__device__ void _warpMin(volatile cudaReal* sData, int tId)
{
   if (sData[tId + 32] < sData[tId]) sData[tId] = sData[tId + 32];
   if (sData[tId + 16] < sData[tId]) sData[tId] = sData[tId + 16];
   if (sData[tId + 8] < sData[tId]) sData[tId] = sData[tId + 8];
   if (sData[tId + 4] < sData[tId]) sData[tId] = sData[tId + 4];
   if (sData[tId + 2] < sData[tId]) sData[tId] = sData[tId + 2];
   if (sData[tId + 1] < sData[tId]) sData[tId] = sData[tId + 1];
}

}
}
