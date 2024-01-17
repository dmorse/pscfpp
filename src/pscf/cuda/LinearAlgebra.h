#ifndef PSCF_LINEAR_ALGEBRA_H
#define PSCF_LINEAR_ALGEBRA_H

#include "GpuTypes.h"
#include <complex>

namespace Pscf {

/** \ingroup Pscf_Cuda_Module 
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

__global__ void inPlacePointwiseDivComplex(cudaComplex* a, const cudaReal* b, int size);

__global__ void pointWiseBinaryMultiply(const cudaReal* a, const cudaReal* b, cudaReal* result, int size);

__global__ void assignUniformReal(cudaReal* result, cudaReal uniform, int size);

__global__ void assignReal(cudaReal* result, const cudaReal* rhs, int size);

__global__ void assignExp(cudaReal* exp, const cudaReal* w, double constant, int size);

__global__ void assignRealSubtractDouble(cudaReal* result, const cudaReal* rhs, double constant, int size);

__global__ void scaleReal(cudaReal* result, double scale, int size);

__global__ void scaleComplex(cudaComplex* result, double scale, int size);

__global__ void mcftsScale(cudaReal* result, cudaReal scale, int size);

__global__ void fourierMove(cudaComplex* a, const cudaReal* b, const cudaReal* c, int size);

__global__ void computeDField(cudaReal* result, const cudaReal* Wc, const cudaReal* Cc, double a, double b, double s, int size);

__global__ void computeForceBias(cudaReal* result, const cudaReal* di, const cudaReal* df, const cudaReal* dwc, double mobility, int size);

/** @} */

}
#endif
