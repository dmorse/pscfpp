#ifndef PSPG_LINEAR_ALGEBRA_H
#define PSPG_LINEAR_ALGEBRA_H

#include "GpuTypes.h"

namespace Pscf {
namespace Pspg {

__global__ void subtractUniform(cudaReal* result, cudaReal rhs, int size);

__global__ void addUniform(cudaReal* result, cudaReal rhs, int size);

__global__ void pointWiseSubtract(cudaReal* result, const cudaReal* rhs, int size);

__global__ void pointWiseSubtractFloat(cudaReal* result, const float rhs, int size);

__global__ void pointWiseBinarySubtract(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void pointWiseBinaryAdd(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void pointWiseBinaryMultiply(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void pointWiseAddScale(cudaReal* result, const cudaReal* rhs, double scale, int size);

__global__ void assignUniformReal(cudaReal* result, cudaReal uniform, int size);

__global__ void assignReal(cudaReal* result, const cudaReal* rhs, int size);

__global__ void inPlacePointwiseMul(cudaReal* a, const cudaReal* b, int size);

}
}
#endif
