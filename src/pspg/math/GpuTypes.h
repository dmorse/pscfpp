#ifndef PSPG_GPU_TYPES_H
#define PSPG_GPU_TYPES_H

#include <cufft.h>
#include <stdio.h>
#include <iostream>

namespace Pscf {
namespace Pspg {

extern int THREADS_PER_BLOCK;
extern int NUMBER_OF_BLOCKS;
extern int MAX_THREADS_PER_BLOCK;

//#define SINGLE_PRECISION
#define DOUBLE_PRECISION


#ifdef SINGLE_PRECISION
typedef cufftReal cudaReal;
typedef cufftComplex cudaComplex;
typedef cufftReal hostReal;
typedef cufftComplex hostComplex;
#else
#ifdef DOUBLE_PRECISION
typedef cufftDoubleReal cudaReal;
typedef cufftDoubleComplex cudaComplex;
typedef cufftDoubleReal hostReal;
typedef cufftDoubleComplex hostComplex;
#endif
#endif

void setGpuBlocksThreads(int & nBlocks, int & nThreads, int const & datasize);

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if(abort) exit(code);
    }
}

}
}
#endif