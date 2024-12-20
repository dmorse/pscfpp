#ifndef PSCF_REDUCE_H
#define PSCF_REDUCE_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include "DeviceArray.h"

namespace Pscf {
namespace Reduce {

/** 
* Functions that perform parallel reductions on the GPU. 
* 
* A reduction is any operation that involves reducing all of the 
* elements of an array (or a set of multiple arrays) down to a single 
* result. Examples include taking the sum or finding the maximum of 
* all array elements. Highly efficient algorithms have been developed 
* to perform such operations in parallel on a GPU, and those 
* algorithms are implemented here.
* 
* A kernel wrapper is provided for each reduction operation, which
* takes a DeviceArray as an input and returns a single value. The 
* wrappers are called on the host, and they call CUDA kernels 
* internally to perform the reductions in parallel.
* 
* The CUDA kernels that perform the parallel reductions are not 
* intended to be called directly. The wrappers and CUDA kernels are 
* implemented together in a way that optimizes the speed of the 
* overall operation. To do this, several important aspects of the
* algorithm are performed by the wrapper, detailed further below.
* 
* First, the CUDA kernels do not perform the entire reduction. 
* Rather, each block of threads in the GPU is reduced to a single 
* value, and the kernel returns an array containing one value for 
* each thread block. The kernel wrapper then performs the remaining 
* reduction by either calling the kernel again or performing the
* remaining reduction on the host, depending on the size of the 
* array output by the first kernel. Once the array is smaller than
* 1e5 elements, the remaining reduction is performed on the CPU. 
* (This value is chosen because the GPU does not significantly speed 
* up reductions on arrays smaller than 1e5 elements.)
* 
* Second, the CUDA kernels are designed to use a number of threads
* that is only half the size of the input array (rounded up to the 
* nearest multiple of the thread block size). The kernel wrapper 
* properly determines the GPU configuration before calling the kernel.
* 
* \ingroup Pscf_Cuda_Module 
* @{
*/

/**
* Compute sum of array elements (GPU kernel wrapper).
*
* \param in  input array
*/
__host__ cudaReal sum(DeviceArray<cudaReal> const & in);

/**
* Get maximum of array elements (GPU kernel wrapper).
*
* \param in  input array
*/
__host__ cudaReal max(DeviceArray<cudaReal> const & in);

/**
* Get maximum absolute magnitude of array elements (GPU kernel wrapper).
*
* \param in  input array
*/
__host__ cudaReal maxAbs(DeviceArray<cudaReal> const & in);

/**
* Get minimum of array elements (GPU kernel wrapper).
*
* \param in  input array
*/
__host__ cudaReal min(DeviceArray<cudaReal> const & in);

/**
* Get minimum absolute magnitude of array elements (GPU kernel wrapper).
*
* \param in  input array
*/
__host__ cudaReal minAbs(DeviceArray<cudaReal> const & in);

/**
* Compute inner product of two real arrays (GPU kernel wrapper).
*
* \param a  first input array
* \param b  second input array
*/
__host__ cudaReal innerProduct(DeviceArray<cudaReal> const & a, 
                               DeviceArray<cudaReal> const & b);

/**
* Compute sum of array elements (GPU kernel).
*
* Assumes each warp is 32 threads. 
* Assumes that each block contains at least 64 threads.
* Assumes that the block size is a power of 2.
*
* \param sum  reduced array containing the sum from each thread block
* \param in  input array
* \param n  number of input array elements
*/
__global__ void _sum(cudaReal* sum, const cudaReal* in, int n);

/**
* Get maximum of array elements (GPU kernel).
*
* Assumes each warp is 32 threads. 
* Assumes that each block contains at least 64 threads.
* Assumes that the block size is a power of 2.
*
* \param max  reduced array containing the max from each thread block
* \param in  input array
* \param n  number of input array elements
*/
__global__ void _max(cudaReal* max, const cudaReal* in, int n);

/**
* Get maximum absolute magnitude of array elements (GPU kernel).
*
* Assumes each warp is 32 threads. 
* Assumes that each block contains at least 64 threads.
* Assumes that the block size is a power of 2.
*
* \param max  reduced array containing the max from each thread block
* \param in  input array
* \param n  number of input array elements
*/
__global__ void _maxAbs(cudaReal* max, const cudaReal* in, int n);

/**
* Get minimum of array elements (GPU kernel).
*
* Assumes each warp is 32 threads. 
* Assumes that each block contains at least 64 threads.
* Assumes that the block size is a power of 2.
*
* \param min  reduced array containing the min from each thread block
* \param in  input array
* \param n  number of input array elements
*/
__global__ void _min(cudaReal* min, const cudaReal* in, int n);

/**
* Get minimum absolute magnitude of array elements (GPU kernel).
*
* Assumes each warp is 32 threads. 
* Assumes that each block contains at least 64 threads.
* Assumes that the block size is a power of 2.
*
* \param min  reduced array containing the min from each thread block
* \param in  input array
* \param n  number of input array elements
*/
__global__ void _minAbs(cudaReal* min, const cudaReal* in, int n);

/**
* Compute inner product of two real arrays (GPU kernel).
*
* Assumes each warp is 32 threads. 
* Assumes that each block contains at least 64 threads.
* Assumes that the block size is a power of 2.
*
* \param ip  reduced array containing the inner prod from each thread block
* \param a  first input array
* \param b  second input array
* \param n  number of input array elements
*/
__global__ void _innerProduct(cudaReal* ip, const cudaReal* a, 
                              const cudaReal* b, int n);

/*
* Parallel reductions performed by a single warp of threads. 
*
* The reduction algorithm can be simplified during the last 6 
* levels of reduction, because these levels of reduction are 
* performed by a single warp. Within a single warp, each thread
* executes the same instruction at the same time (SIMD execution). 
* Therefore, we don't need the __syncthreads() command between 
* reduction operations. Further, we do not need to evaluate an 
* if-statement to determine which threads should perform the 
* calculation and which should not, since the entire warp will be 
* dedicated to these operations regardless of whether they perform 
* calculations. Therefore, the if-statement would not free any
* resources for other tasks, so we omit it for speed.
* 
* We assume here that a single warp contains 32 threads. All 
* CUDA-compatible GPUs currently meet this criterion, but it is 
* possible that someday there will be GPUs with a different warp
* size. The methods below may break if the warp size is smaller
* than 32 threads, because the operations would be performed by
* multiple warps without __syncthreads() commands to keep them
* synced. Warps larger than 32 threads would still be compatible
* with these functions, though the functions are not optimized
* for this case. 
* 
* These are implemented as separate functions, rather than within
* the kernels above, because they require the sData array to be
* defined as volatile (meaning the array values may change at any
* time, so the compiler must access the actual memory location 
* rather than using cached values).
*/

/**
* Utility to perform parallel reduction summation within a single warp.
*
* \param sData  input array to reduce
* \param tId  thread ID
*/
__device__ void _warpSum(volatile cudaReal* sData, int tId);

/**
* Utility to perform parallel reduction maximization within a single warp.
*
* \param sData  input array to reduce
* \param tId  thread ID
*/
__device__ void _warpMax(volatile cudaReal* sData, int tId);

/**
* Utility to perform parallel reduction minimization within a single warp.
*
* \param sData  input array to reduce
* \param tId  thread ID
*/
__device__ void _warpMin(volatile cudaReal* sData, int tId);

/** @} */

}
}
#endif
