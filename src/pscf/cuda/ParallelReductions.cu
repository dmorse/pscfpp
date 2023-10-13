#ifndef PSCF_PARALLEL_REDUCTIONS_CU
#define PSCF_PARALLEL_REDUCTIONS_CU

#include "ParallelReductions.h"

namespace Pscf {

__global__ void reductionSum(cudaReal* sum, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int bdim = blockDim.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Global memory load and first operation
      cudaReal in0 = in[idx];
      if (idx + bdim < size) {
         cudaReal in1 = in[idx+bdim];
         sdata[tid] = in0+in1;
      } else {
         sdata[tid] = in0;
      }
      
      // wait for all threads to finish
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride && idx+stride < size) {
            sdata[tid] += sdata[tid+stride];
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (tid == 0)
         sum[bid] = sdata[0];
   }
}

__global__ void reductionInnerProduct(cudaReal* innerprod, const cudaReal* a, const cudaReal* b, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int bdim = blockDim.x;
   int idx = bid * bdim*2 + tid;

   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Global memory load and first operation
      cudaReal in0 = a[idx]*b[idx];
      if (idx + bdim < size) {
         cudaReal in1 = a[idx+bdim]*b[idx+bdim];
         sdata[tid] = in0 + in1;
      } else {
         sdata[tid] = in0;
      }
      
      
      // wait for all threads to finish
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride && idx+stride < size) {
            sdata[tid] += sdata[tid+stride];
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (tid == 0) {
         innerprod[bid] = sdata[0];
      }
   }
}

__global__ void reductionMax(cudaReal* max, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int bdim = blockDim.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Global memory load and first operation
      cudaReal in0 = in[idx];
      if (idx + bdim < size) {
         cudaReal in1 = in[idx+bdim];
         sdata[tid] = (in0 > in1) ? in0 : in1;
      } else {
         sdata[tid] = in0;
      }
      
      // wait for all threads to finish
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride && idx+stride < size) {
            if (sdata[tid+stride] > sdata[tid]) {
               sdata[tid] = sdata[tid+stride];
            }
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (tid == 0)
         max[bid] = sdata[0];
   }
}

__global__ void reductionMaxAbs(cudaReal* max, const cudaReal* in, int size) 
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int bdim = blockDim.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Global memory load and first operation
      cudaReal in0 = fabs(in[idx]);
      if (idx + bdim < size) {
         cudaReal in1 = fabs(in[idx+bdim]);
         sdata[tid] = (in0 > in1) ? in0 : in1;
      } else {
         sdata[tid] = in0;
      }
      
      // wait for all threads to finish
      __syncthreads();
      
      // reduction
      // data and block dimensions need to be a power of two
      for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
         if (tid < stride && idx+stride < size) {
            if (sdata[tid+stride] > sdata[tid]) {
                  sdata[tid] = sdata[tid+stride];
               }
         }
         __syncthreads();
      }
         
      // one thread for each block stores that block's results in global memory
      if (tid == 0)
         max[bid] = sdata[0];
   }
}

__global__ void reductionMin(cudaReal* min, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int bdim = blockDim.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Global memory load and first operation
      cudaReal in0 = in[idx];
      if (idx + bdim < size) {
         cudaReal in1 = in[idx+bdim];
         sdata[tid] = (in0 < in1) ? in0 : in1;
      } else {
         sdata[tid] = in0;
      }
      
      // wait for all threads to finish
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride && idx+stride < size) {
            if (sdata[tid+stride] < sdata[tid]) {
               sdata[tid] = sdata[tid+stride];
            }
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (tid == 0) {
         min[bid] = sdata[0];
      }
   }
}

__global__ void reductionMinAbs(cudaReal* min, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int bdim = blockDim.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Global memory load and first operation
      cudaReal in0 = fabs(in[idx]);
      if (idx + bdim < size) {
         cudaReal in1 = fabs(in[idx+bdim]);
         sdata[tid] = (in0 < in1) ? in0 : in1;
      } else {
         sdata[tid] = in0;
      }
      
      // wait for all threads to finish
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride && idx+stride < size) {
            if (sdata[tid+stride] < sdata[tid]) {
               sdata[tid] = sdata[tid+stride];
            }
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (tid == 0) {
         min[bid] = sdata[0];
      }
   }
}

}
#endif
