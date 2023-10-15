#ifndef PSCF_KERNEL_WRAPPERS_H
#define PSCF_KERNEL_WRAPPERS_H

#include "ParallelReductions.h"
#include "GpuTypes.h"
#include "LinearAlgebra.h"

namespace Pscf {

/** \ingroup Pscf_Cuda_Module 
* @{
*/

/**
* Compute sum of array elements (CUDA host function).
*
* \param d_in  input array
* \param size  size of array (number of elements).
*/
__host__ cudaReal gpuSum(cudaReal const * d_in, int size);

/**
* Compute maximum of array elements (CUDA host function).
*
* \param d_in  input array
* \param size  size of array (number of elements).
*/
__host__ cudaReal gpuMax(cudaReal const * d_in, int size);

/**
* Compute maximum absolute mag. of array elements (CUDA host function).
*
* \param d_in  input array
* \param size  size of array (number of elements).
*/
__host__ cudaReal gpuMaxAbs(cudaReal const * d_in, int size);

/**
* Compute minimum of array elements (CUDA host function).
*
* \param d_in  input array
* \param size  size of array (number of elements)
*/
__host__ cudaReal gpuMin(cudaReal const * d_in, int size);

/**
* Compute minimum absolut mag. of array elements (CUDA host function).
*
* \param d_in  input array
* \param size  size of array (number of elements)
*/
__host__ cudaReal gpuMinAbs(cudaReal const * d_in, int size);

/**
* Compute inner product of two real arrays (CUDA host function).
*
* \param d_a  first input array
* \param d_b  second input array
* \param size  size of array (number of elements)
*/
__host__ cudaReal 
gpuInnerProduct(cudaReal const * d_a, cudaReal const * d_b, int size);

/** @} */

}
#endif
