#ifndef RPG_HELPERS_CU
#define RPG_HELPERS_CU

#include "Helpers.h"

namespace Pscf {
namespace Rpg {

__global__ void helmholtzHelper(cudaReal* result, const cudaReal* composition,
   const cudaReal* pressure, double chi, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = composition[i] * composition[i] / chi - pressure[i];
   }
}

__global__ void reformField(cudaReal* Wa, cudaReal* Wb,
   const cudaReal* pressureF, const cudaReal* compositionF, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      Wa[i] = pressureF[i] + compositionF[i];
      Wb[i] = pressureF[i] - compositionF[i];
   }
}

__global__ void mcftsStepHelper(cudaReal* result, const cudaReal* A2, const cudaReal* sqrtSq, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
          result[i] += A2[i] * sqrtSq[i];
   }
}
__global__ void mcftsScale(cudaReal* result, cudaReal scale, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
          result[i] = result[i] * 2 * scale - scale;
   }
}

}
}
#endif
