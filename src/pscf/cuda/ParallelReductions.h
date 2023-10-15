#ifndef PARALLEL_REDUCTIONS_H
#define PARALLEL_REDUCTIONS_H

#include "GpuTypes.h"

namespace Pscf {

/** \ingroup Pscf_Cuda_Module 
* @{
*/
__global__ void reductionSum(cudaReal* sum, const cudaReal* in, int size);

__global__ void reductionInnerProduct(cudaReal* innerprod, const cudaReal* a, const cudaReal* b, int size);

__global__ void reductionMax(cudaReal* max, const cudaReal* in, int size);

__global__ void reductionMaxAbs(cudaReal* max, const cudaReal* in, int size);

__global__ void reductionMin(cudaReal* min, const cudaReal* in, int size);

__global__ void reductionMinAbs(cudaReal* min, const cudaReal* in, int size);

/** @} */

}
#endif
