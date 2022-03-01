#ifndef PSPG_KERNEL_WRAPPERS_H
#define PSPG_KERNEL_WRAPPERS_H

#include "ParallelReductions.h"
#include "GpuTypes.h"
#include "LinearAlgebra.h"

namespace Pscf {
namespace Pspg {

/** \ingroup Pspg_Math_Module 
* @{
*/
__host__ cudaReal gpuSum(cudaReal const * d_in, int size);

__host__ cudaReal gpuInnerProduct(cudaReal const * d_a, cudaReal const * d_b, int size);

__host__ cudaReal gpuMax(cudaReal const * d_in, int size);

__host__ cudaReal gpuMaxAbs(cudaReal const * d_in, int size);

__host__ cudaReal gpuMin(cudaReal const * d_in, int size);

__host__ cudaReal gpuMinAbs(cudaReal const * d_in, int size);

/** @} */

}
}
#endif
