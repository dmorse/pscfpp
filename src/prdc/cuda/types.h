#ifndef PRDC_CUDA_TYPES_H
#define PRDC_CUDA_TYPES_H

#include <cufft.h>
#include <stdio.h>
#include <iostream>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   /**
   * Complex number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   typedef cufftComplex Complex;
   #else
   #ifdef DOUBLE_PRECISION
   typedef cufftDoubleComplex Complex;
   #endif
   #endif

   /**
   * Real number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   typedef cufftReal Real;
   #else
   #ifdef DOUBLE_PRECISION
   typedef cufftDoubleReal Real;
   #endif
   #endif

} // Pscf::Prdc::Cpu
} // Pscf::Prdc
} // Pscf
#endif
