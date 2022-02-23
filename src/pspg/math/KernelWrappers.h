#ifndef PSPG_KERNEL_WRAPPERS_H
#define PSPG_KERNEL_WRAPPERS_H

#include "ParallelReductions.h"
#include "GpuTypes.h"
#include "LinearAlgebra.h"

namespace Pscf {
namespace Pspg {

// __host__ cudaReal gpuSum(cudaReal* d_in, int size);

__host__ cudaReal gpuInnerProduct(cudaReal* d_a, cudaReal* d_b, int size);

// __host__ cudaReal gpuMax(cudaReal* d_in, int size);

// __host__ cudaReal gpuMaxAbs(cudaReal* d_in, int size);

// __host__ cudaReal gpuMin(cudaReal* d_in, int size);

// __host__ cudaReal gpuMinAbs(cudaReal* d_in, int size);

}
}
#endif
