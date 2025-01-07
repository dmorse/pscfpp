#ifndef PSCF_GPU_TYPES_H
#define PSCF_GPU_TYPES_H

#include <cufft.h>
#include <stdio.h>
#include <iostream>

namespace Pscf {

//#define SINGLE_PRECISION
#define DOUBLE_PRECISION

#ifdef SINGLE_PRECISION
typedef cufftReal cudaReal;
typedef cufftComplex cudaComplex;
#else
#ifdef DOUBLE_PRECISION
typedef cufftDoubleReal cudaReal;
typedef cufftDoubleComplex cudaComplex;
#endif
#endif

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, 
                      bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", 
              cudaGetErrorString(code), file, line);
      if(abort) exit(code);
    }
}

}
#endif
