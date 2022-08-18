#ifndef PSPG_HELPERS_H
#define PSPG_HELPERS_H

#include "GpuTypes.h"

namespace Pscf {
namespace Pspg {

__global__ void helmholtzHelper(cudaReal* result, const cudaReal* composition,
   const cudaReal* pressure, double chi, int size);

__global__ void reformField(cudaReal* Wa, cudaReal* Wb,
   const cudaReal* pressureF,const cudaReal* compositionF, int size);

__global__ void mcftsStepHelper(cudaReal* result, const cudaReal* A2, const cudaReal* sqrtSq, int size);

__global__ void mcftsScale(cudaReal* result, cudaReal scale, int size);

}
}
#endif
