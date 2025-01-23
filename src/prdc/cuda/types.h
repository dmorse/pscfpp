#ifndef PRDC_CUDA_TYPES_H
#define PRDC_CUDA_TYPES_H

#include <cufft.h>

namespace Pscf {

   // Toggle single / double precision:
   
   //#define SINGLE_PRECISION
   #define DOUBLE_PRECISION

   /**
   * Complex number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   typedef cufftComplex Complex;
   typedef cufftComplex cudaComplex;
   #else
   #ifdef DOUBLE_PRECISION
   typedef cufftDoubleComplex Complex;
   typedef cufftDoubleComplex cudaComplex;
   #endif
   #endif

   /**
   * Real number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   typedef cufftReal Real;
   typedef cufftReal cudaReal;
   #else
   #ifdef DOUBLE_PRECISION
   typedef cufftDoubleReal Real;
   typedef cufftDoubleReal cudaReal;
   #endif
   #endif

} // Pscf
#endif
