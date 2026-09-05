#ifndef PSCF_CUDA_TYPES_H
#define PSCF_CUDA_TYPES_H

#include <cufft.h>
//#include <cuComplex.h>

// Toggle single / double precision:
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//#define SINGLE_PRECISION
#define DOUBLE_PRECISION

namespace Pscf {

   /**
   * Complex number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   using cudaComplex = cufftComplex;
   //using cudaComplex = cuComplex;
   #else
   #ifdef DOUBLE_PRECISION
   using cudaComplex = cufftDoubleComplex;
   //using cudaComplex = cuDoubleComplex;
   #endif
   #endif

   /**
   * Real number type used in CPU code that uses FFTW.
   */
   #ifdef SINGLE_PRECISION
   using cudaReal = cufftReal;
   //using cudaReal = float;
   #else
   #ifdef DOUBLE_PRECISION
   using cudaReal = cufftDoubleReal;
   //using cudaReal = double;
   #endif
   #endif

}
#endif
