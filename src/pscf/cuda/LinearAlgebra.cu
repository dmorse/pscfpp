#include "LinearAlgebra.h"

namespace Pscf 
{

__global__ void subtractUniform(cudaReal* result, cudaReal rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs;
   }
}

__global__ void addUniform(cudaReal* result, cudaReal rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += rhs;
   }
}

__global__ void pointWiseSubtract(cudaReal* result, const cudaReal* rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs[i];
   }
}

__global__ void pointWiseSubtractFloat(cudaReal* result, const float rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs;
   }   
}

__global__ void pointWiseBinarySubtract(const cudaReal* a, const cudaReal* b, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] - b[i];
   }
}

__global__ void pointWiseAdd(cudaReal* result, const cudaReal* rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += rhs[i];
   }
}

__global__ void pointWiseBinaryAdd(const cudaReal* a, const cudaReal* b, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] + b[i];
   }
}

__global__ void pointWiseAddScale(cudaReal* result, const cudaReal* rhs, double scale, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += scale * rhs[i];
   }
}

__global__ void inPlacePointwiseMul(cudaReal* a, const cudaReal* b, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      a[i] *= b[i];
   }
}

__global__ void pointWiseBinaryMultiply(const cudaReal* a, const cudaReal* b, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] * b[i];
   }
}

__global__ void assignUniformReal(cudaReal* result, cudaReal uniform, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = uniform;
   }
}

__global__ void assignReal(cudaReal* result, const cudaReal* rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = rhs[i];
   }
}

__global__ void assignExp(cudaReal* out, const cudaReal* w, double constant, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      out[i] = exp(-w[i]*constant);
   }
}

__global__ void scaleReal(cudaReal* result, double scale, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;

   for (int i = startID; i < size; i += nThreads) {
      result[i] *= scale;
   }
}

__global__ void mcftsScale(cudaReal* result, cudaReal scale, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      // Scale result from [0,1] to [-scale, scale]
      result[i] = result[i] * 2 * scale - scale;
   }
}

__global__ void fourierMove(cudaComplex* a, const cudaReal* b, const cudaReal* c, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      a[i].x += b[i];
      a[i].y += c[i];
   }
}

}
