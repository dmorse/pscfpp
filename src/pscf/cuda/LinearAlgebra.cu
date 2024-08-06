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

__global__ void inPlacePointwiseDivComplex(cudaComplex* a, const cudaReal* b, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      a[i].x /= b[i];
      a[i].y /= b[i];
   }
}

__global__ void squaredMagnitudeComplex(const cudaComplex* a, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i].x *  a[i].x + a[i].y * a[i].y;
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

__global__ void assignRealSubtractDouble(cudaReal* result, const cudaReal* rhs, double constant, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = rhs[i] - constant;
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

__global__ void scaleComplex(cudaComplex* result, double scale, int size){
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i].x *= scale;
      result[i].y *= scale;
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

__global__ void computeDField(cudaReal* result, const cudaReal* Wc, const cudaReal* Cc, double a, double b, double s, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a*( b*(Wc[i] - s) + Cc[i] );
   }
}

__global__ void computeForceBias(cudaReal* result, const cudaReal* di, const cudaReal* df, const cudaReal* dwc, double mobility, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = 0.5*(di[i] + df[i]) * (dwc[i] + mobility* (0.5*(di[i] - df[i])));
   }
}


}
