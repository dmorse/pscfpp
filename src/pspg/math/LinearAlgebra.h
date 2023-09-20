#ifndef PSPG_LINEAR_ALGEBRA_H
#define PSPG_LINEAR_ALGEBRA_H

#include "GpuTypes.h"
#include <complex>

namespace Pscf {
namespace Pspg {

/** \ingroup Pspg_Math_Module 
* @{
*/


__global__ void subtractUniform(cudaReal* result, cudaReal rhs, int size);

__global__ void addUniform(cudaReal* result, cudaReal rhs, int size);

__global__ void pointWiseSubtract(cudaReal* result, const cudaReal* rhs, int size);

__global__ void pointWiseSubtractFloat(cudaReal* result, const float rhs, int size);

__global__ void pointWiseBinarySubtract(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void pointWiseAdd(cudaReal* result, const cudaReal* rhs, int size);

__global__ void pointWiseBinaryAdd(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void pointWiseAddScale(cudaReal* result, const cudaReal* rhs, double scale, int size);

__global__ void inPlacePointwiseMul(cudaReal* a, const cudaReal* b, int size);

__global__ void pointWiseBinaryMultiply(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void assignUniformReal(cudaReal* result, cudaReal uniform, int size);

__global__ void assignReal(cudaReal* result, const cudaReal* rhs, int size);

__global__ void assignExp(cudaReal* exp, const cudaReal* w, double constant, int size);

__global__ void scaleReal(cudaReal* result, double scale, int size);

__global__ void mcftsScale(cudaReal* result, cudaReal scale, int size);

__global__ void fourierMove(cudaComplex* a, const cudaReal* b, const cudaReal* c, int size);

/** @} */

}
}
#endif
