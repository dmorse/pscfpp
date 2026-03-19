/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOp.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/cudaErrorCheck.h>
#include <cmath>

namespace Pscf {
namespace VecOp {

// CUDA kernels:
// (defined in anonymous namespace, used within this file only)

namespace {

      #if 0
      /*
      * Vector assignment, a[i] = b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i];
         }
      }
      #endif

      /*
      * Vector assignment, a[i] = b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x;
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector assignment, a[j] = (b[j], c[j]) (complex, real & imaginary).
      *
      * \param a  complex array (LHS)
      * \param b  real array, real part (RHS)
      * \param c  real array, imaginary part (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqV(cudaComplex* a, 
                cudaReal const * b, cudaReal const * c, 
                const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i];
            a[i].y = c[i];
         }
      }

      /*
      * Vector assignment, a[i] = b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b;
         }
      }

      /*
      * Vector assignment, a[i] = b(complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqS(cudaComplex* a, const cudaComplex b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b.x;
            a[i].y = b.y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaReal* a,
                  cudaReal const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] + c[i];
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i](complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaComplex* a,
                  cudaComplex const * b,
                  cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c[i].x;
            a[i].y = b[i].y + c[i].y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i] (mixed, b real).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaComplex* a,
                  cudaReal const * b,
                  cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] + c[i].x;
            a[i].y = c[i].y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c[i] (mixed, c real).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVV(cudaComplex* a,
                  cudaComplex const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c[i];
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaReal* a,
                  cudaReal const * b,
                  const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] + c;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c(complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaComplex* a,
                  cudaComplex const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c.x;
            a[i].y = b[i].y + c.y;
         }
      }

      /*
      * Vector addition, a[i] = b[i] + c (mixed, b real).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaComplex* a,
                  cudaReal const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] + c.x;
            a[i].y = c.y;
         }
      }

      /*
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVS(cudaComplex* a,
                  cudaComplex const * b,
                  const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x + c;
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param c  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subVV(cudaReal* a,
                  cudaReal const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] - c[i];
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i](complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVV(cudaComplex* a, cudaComplex const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c[i].x;
            a[i].y = b[i].y - c[i].y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i], (mixed, b real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _subVV(cudaComplex* a, cudaReal const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] - c[i].x;
            a[i].y = 0.0 - c[i].y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c[i] (mixed, c real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _subVV(cudaComplex* a, 
                  cudaComplex const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c[i];
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector-scalar subtraction, a[i] = b[i] - c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _subVS(cudaReal* a, 
                  cudaReal const * b,
                  const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] - c;
         }
      }

      /*
      * Vector-scalar subtraction, a[i] = b[i] - c (complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _subVS(cudaComplex* a, 
                  cudaComplex const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c.x;
            a[i].y = b[i].y - c.y;
         }
      }

      /*
      * Vector-scalar subtraction, a[i] = b[i] - c (mixed, b real).
      *
      * \param a  complex output array (LHS)
      * \param b  real input array (RHS)
      * \param c  real input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _subVS(cudaComplex* a, cudaReal const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] - c.x;
            a[i].y = 0.0 - c.y;
         }
      }

      /*
      * Vector subtraction, a[i] = b[i] - c (mixed, c real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _subVS(cudaComplex* a, 
                  cudaComplex const * b,
                  const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x - c;
            a[i].y = b[i].y;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _mulVV(cudaReal* a, 
                  cudaReal const * b,
                  cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] * c[i];
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i](complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _mulVV(cudaComplex* a, 
                  cudaComplex const * b,
                  cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = (b[i].x * c[i].x) - (b[i].y * c[i].y);
            a[i].y = (b[i].x * c[i].y) + (b[i].y * c[i].x);
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i] (mixed, b real).
      *
      * \param a  complex output array (LHS)
      * \param b  real input array (RHS)
      * \param c  complex input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVV(cudaComplex* a, 
                             cudaReal const * b,
                             cudaComplex const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] * c[i].x;
            a[i].y = b[i] * c[i].y;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c[i] (mixed, c real).
      *
      * \param a  complex output array (LHS)
      * \param b  complex input array (RHS)
      * \param c  real input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVV(cudaComplex* a, 
                             cudaComplex const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x * c[i];
            a[i].y = b[i].y * c[i];
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVS(cudaReal* a, cudaReal const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] * c;
         }
      }

      /*
      * Vector-scalar multiplication, a[i] = b[i] * c(complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _mulVS(cudaComplex* a, 
                  cudaComplex const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = (b[i].x * c.x) - (b[i].y * c.y);
            a[i].y = (b[i].x * c.y) + (b[i].y * c.x);
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c (mixed, b real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _mulVS(cudaComplex* a, 
                  cudaReal const * b,
                  const cudaComplex c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i] * c.x;
            a[i].y = b[i] * c.y;
         }
      }

      /*
      * Vector multiplication, a[i] = b[i] * c (mixed, c real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _mulVS(cudaComplex* a, cudaComplex const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x * c;
            a[i].y = b[i].y * c;
         }
      }

      /*
      * Vector division, a[i] = b[i] / c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVV(cudaReal* a, cudaReal const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] / c[i];
         }
      }

      /*
      * Vector division, a[i] = b[i] / c[i] (mixed, c real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVV(cudaComplex* a, cudaComplex const * b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x / c[i];
            a[i].y = b[i].y / c[i];
         }
      }

      /*
      * Vector division, a[i] = b[i] / c (real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVS(cudaReal* a, cudaReal const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] / c;
         }
      }

      /*
      * Vector division, a[i] = b[i] / c (mixed, c real).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__ void _divVS(cudaComplex* a, cudaComplex const * b,
                             const cudaReal c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = b[i].x / c;
            a[i].y = b[i].y / c;
         }
      }

      /*
      * Vector division, a[i] = b / c[i] (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param c  input array (RHS)
      * \param n  size of arrays
      */
      __global__ void _divSV(cudaReal* a, const cudaReal b,
                             cudaReal const * c, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b / c[i];
         }
      }

      // In-place addition

      /*
      * Vector in-place addition, a[i] += b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b[i];
         }
      }

      /*
      * Vector in-place addition, a[i] += b[i] (complex).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _addEqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i].x;
            a[i].y += b[i].y;
         }
      }

      /*
      * Vector in-place addition, a[i] += (b[i],c[i]) (complex, real/imag).
      *
      * \param a  complex input / output array (LHS)
      * \param b  real array, increment to real part (RHS)
      * \param c  real array, increment to imaginary part (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _addEqV(cudaComplex* a,
                   cudaReal const * b, 
                   cudaReal const * c, 
		   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i];
            a[i].y += c[i];
         }
      }

      /*
      * Vector addition in-place, a[i] += b[i] (mixed).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param n  size of arrays
      */
      __global__ 
      void _addEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b[i];
         }
      }

      /*
      * Vector addition in-place, a[i] += b (real).
      *
      * \param a  output array (LHS)
      * \param b  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b;
         }
      }

      /*
      * Vector-scalar in-place addition, a[i] += b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqS(cudaComplex* a, const cudaComplex b,
                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b.x;
            a[i].y += b.y;
         }
      }

      /*
      * Vector addition in-place, a[i] += b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x += b;
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] -= b[i];
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b[i] (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b[i].x;
            a[i].y -= b[i].y;
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b[i];
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] -= b;
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqS(cudaComplex* a, const cudaComplex b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b.x;
            a[i].y -= b.y;
         }
      }

      /*
      * Vector in-place subtraction, a[i] -= b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x -= b;
         }
      }

      /*
      * Vector in-place multiplication, a[i] *= b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] *= b[i];
         }
      }

      /*
      * Vector in-place multiplication, a[i] *= b[i], (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaComplex c;
         for (int i = startID; i < n; i += nThreads) {
            c.x = (a[i].x * b[i].x) - (a[i].y * b[i].y);
            c.y = (a[i].x * b[i].y) + (a[i].y * b[i].x);
            a[i].x = c.x;
            a[i].y = c.y;
         }
      }

      /*
      * Vector in-place multiplication, a[i]*=b[i] (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x *= b[i];
            a[i].y *= b[i];
         }
      }

      /*
      * Vector-scalar in-place multiplication, a[i] *= b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] *= b;
         }
      }

      /*
      * Vector-scalar in-place multiplication, a[i] *= b (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqS(cudaComplex* a, const cudaComplex b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaComplex c;
         for (int i = startID; i < n; i += nThreads) {
            c.x = (a[i].x * b.x) - (a[i].y * b.y);
            c.y = (a[i].x * b.y) + (a[i].y * b.x);
            a[i].x = c.x;
            a[i].y = c.y;
         }
      }

      /*
      * Vector in-place multiplication, a[i] *= b (mixed).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _mulEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x *= b;
            a[i].y *= b;
         }
      }

      /*
      * Vector in-place division, a[i] /= b[i] (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] /= b[i];
         }
      }

      /*
      * Vector in-place division, a[i] /= b[i] (mixed, b = real).
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqV(cudaComplex* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x /= b[i];
            a[i].y /= b[i];
         }
      }

      /*
      * Vector in-place division, a[i] /= b (real).
      *
      * \param a  real array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqS(cudaReal* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] /= b;
         }
      }

      /*
      * Vector in-place division, a[i] /= b (mixed, b = real).
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqS(cudaComplex* a, const cudaReal b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x /= b;
            a[i].y /= b;
         }
      }

      // Exponentiation

      /*
      * Vector exponentiation, a[i] = exp(b[i]) (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _expV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = exp(b[i]);
         }
      }

      /*
      * Vector exponentiation, a[i] = exp(b[i]) (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _expV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x = exp(b[i].x) * cos(b[i].y);
            a[i].y = exp(b[i].x) * sin(b[i].y);
         }
      }

      // Square

      /*
      * Square of real vector, a[i] = (b[i])^2 (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _sqV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i]*b[i];
         }
      }

      /*
      * Square of complex vector, a[i] = (b[i])^2 (complex).
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _sqV(cudaComplex* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal bx, by;
         for (int i = startID; i < n; i += nThreads) {
            bx = b[i].x;
            by = b[i].y;
            a[i].x = (bx*bx) - (by*by);
            a[i].y = 2.0 * bx * by;
         }
      }

      // Absolute magnitude

      /*
      * Vector absolute magnitude, a[i] = abs(b[i]) (real).
      *
      * \param a  real array (LHS)
      * \param b  real array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _absV(cudaReal* a, cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = std::fabs(b[i]);
         }
      }

      /*
      * Vector absolute magnitude, a[i] = |b[i]|^2 (complex).
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _sqAbsV(cudaReal* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal bx, by;
         for (int i = startID; i < n; i += nThreads) {
            bx = b[i].x;
            by = b[i].y;
            a[i] = bx*bx + by*by;
         }
      }

   } // end anonymous namespace

   // CUDA kernel wrappers:

   /*
   * Vector assignment, a[i] = b[i] (real).
   */
   void eqV(DeviceArray<cudaReal>& a,
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      #if 0
      // Launch kernel
      _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                  b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
      #endif

      // Copy all elements
      cudaErrorCheck( cudaMemcpy(a.cArray() + beginIdA, 
			         b.cArray() + beginIdB, 
                                 n * sizeof(cudaReal), 
                                 cudaMemcpyDeviceToDevice) );

   }

   /*
   * Vector assignment, a[i] = b[i] (real, device to host).
   */
   void eqV(Array<cudaReal>& a,
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Copy all elements
      cudaErrorCheck( cudaMemcpy(a.cArray() + beginIdA, 
			         b.cArray() + beginIdB, 
                                 n * sizeof(cudaReal), 
                                 cudaMemcpyDeviceToHost) );

   }

   /*
   * Vector assignment, a[i] = b[i] (real, device to host).
   */
   void eqV(DeviceArray<cudaReal>& a,
            Array<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Copy all elements
      cudaErrorCheck( cudaMemcpy(a.cArray() + beginIdA, 
			         b.cArray() + beginIdB, 
                                 n * sizeof(cudaReal), 
                                 cudaMemcpyHostToDevice) );

   }

   /*
   * Vector assignment, a[i] = b[i] (complex).
   */
   void eqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                  b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector assignment, a[i] = (b[i], c[i]) (complex, real & imaginary).
   */
   void eqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaReal> const & b,  // real part
            DeviceArray<cudaReal> const & c,  // imaginary part
            const int beginIdA, const int beginIdB, const int beginIdC,
	    const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                  b.cArray()+beginIdB, 
                                  c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar assignment, a[i] = b (real).
   */
   void eqS(DeviceArray<cudaReal>& a,
            const cudaReal b,
            const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar assignment, a[i] = b (complex).
   */
   void eqS(DeviceArray<cudaComplex>& a,
            const cudaComplex b,
            const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Binary arithmetic operations, separate array for result

   /*
   * Vector addition, a[i] = b[i] + c[i] (real).
   */
   void addVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector addition, a[i] = b[i] + c[i] (complex).
   */
   void addVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector addition, a[i] = b[i] + c[i] (mixed, b real).
   */
   void addVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, 
                                    b.cArray() + beginIdB,
                                    c.cArray() + beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector addition, a[i] = b[i] + c[i] (mixed, c real).
   */
   void addVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, 
                                    b.cArray() + beginIdB,
                                    c.cArray() + beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar addition, a[i] = b[i] + c (real).
   */
   void addVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar addition, a[i] = b[i] + c (complex).
   void addVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition, a[i] = b[i] + c (mixed, b real).
   void addVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector addition, a[i] = b[i] + c (mixed, c real).
   void addVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c[i] (real).
   void subVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA,
              const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c[i] (cudaComplex).
   void subVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i]=b[i]-c[i] (mixed).
   void subVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i]=b[i]-c[i] (mixed).
   void subVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c (cudaReal).
   void subVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
               const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector subtraction, a[i] = b[i] - c (complex).
   void subVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaComplex c, const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar subtraction, a[i] = b[i] - c (mixed).
   void subVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              const cudaComplex c, const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar subtraction, a[i] = b[i] - c (mixed).
   void subVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i] = b[i] * c[i] (real).
   void mulVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i] = b[i] * c[i] (cudaComplex).
   void mulVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i]=b[i]*c[i] (mixed).
   void mulVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaComplex> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i]=b[i]*c[i] (mixed).
   void mulVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar multiplication, a[i] = b[i] * c (real).
   void mulVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar multiplication, a[i] = b[i] * c (complex).
   void mulVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector-scalar multiplication, a[i] = b[i] * c (mixed).
   void mulVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaReal> const & b,
              const cudaComplex c,
              const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector multiplication, a[i] = b[i] * c (mixed).
   void mulVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector division, a[i] = b[i] / c[i] (real).
   void divVV(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              DeviceArray<cudaReal> const & c, const int beginIdA,
              const int beginIdB, const int beginIdC, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB,
                                    c.cArray()+beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector division, a[i] = b[i] / c[i] (mixed).
   */
   void divVV(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdB, const int beginIdC, 
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, 
                                    b.cArray() + beginIdB,
                                    c.cArray() + beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar division, a[i] = b[i] / c (real).
   */
   void divVS(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b,
              const cudaReal c,
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, 
                                    b.cArray() + beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector division, a[i] = b[i] / c (mixed, c real).
   */
   void divVS(DeviceArray<cudaComplex>& a,
              DeviceArray<cudaComplex> const & b,
              const cudaReal c, 
              const int beginIdA, const int beginIdB,
              const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                    b.cArray()+beginIdB,
                                    c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Division of scalar by  vector, a[i] = b / c[i] (real).
   */
   void divSV(DeviceArray<cudaReal>& a,
              const cudaReal b,
              DeviceArray<cudaReal> const & c,
              const int beginIdA, const int beginIdC, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divSV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b,
                                    c.cArray() + beginIdC, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // In-place arithmetic operations

   /*
   * Vector addition in-place, a[i] += b[i] (real).
   */
   void addEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector addition in-place, a[i] += b[i] (complex).
   */
   void addEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector addition in-place, a[i] += (b[i], c[i]) (complex, real/imag).
   */
   void addEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               DeviceArray<cudaReal> const & c,
               const int beginIdA, const int beginIdB, const int beginIdC,
	       const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);
      UTIL_CHECK(c.capacity() >= n + beginIdC);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, 
                                     c.cArray() + beginIdB, 
				     n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector addition in-place, a[i] += b[i] (mixed).
   */
   void addEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, 
	       const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place addition, a[i] += b (real).
   */
   void addEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place addition, a[i] += b (complex).
   */
   void addEqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place addition, a[i] += b (mixed).
   */
   void addEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector in-place subtraction, a[i] -= b[i] (real).
   */
   void subEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector in-place subtraction, a[i] -= b[i] (complex).
   */
   void subEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector in-place subtraction, a[i] -= b[i] (mixed).
   */
   void subEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place subtraction, a[i] -= b (real).
   */
   void subEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place subtraction, a[i] -= b (complex).
   */
   void subEqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place subtraction, a[i] -= b (mixed).
   */
   void subEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _subEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // In-place multiplication

   // Vector in-place multiplication, a[i] *= b[i] (real).
   void mulEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector in-place multiplication, a[i] *= b[i] (cudaComplex).
   void mulEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA,
                                     b.cArray()+beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector in-place multiplication, a[i] *= b[i] (mixed).
   void mulEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector in-place multiplication, a[i] *= b (cudaReal).
   void mulEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Vector in-place multiplication, a[i] *= b (cudaComplex).
   void mulEqS(DeviceArray<cudaComplex>& a,
               const cudaComplex b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector in-place multiplication, a[i] *= b (mixed).
   */
   void mulEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // In-place division

   /*
   * Vector elementwise in-place division, a[i] /= b[i] (real).
   */
   void divEqV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector elementwise in-place division, a[i] /= b[i] (mixed).
   */
   void divEqV(DeviceArray<cudaComplex>& a,
               DeviceArray<cudaReal> const & b,
               const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place division, a[i] /= b (real).
   */
   void divEqS(DeviceArray<cudaReal>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector-scalar in-place division, a[i] /= b (mixed).
   */
   void divEqS(DeviceArray<cudaComplex>& a,
               const cudaReal b,
               const int beginIdA, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqS<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Exponentiation

   /*
   * Vector exponentiation, a[i] = exp(b[i]) (real).
   */
   void expV(DeviceArray<cudaReal>& a,
             DeviceArray<cudaReal> const & b,
             const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _expV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                   b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector exponentiation, a[i] = exp(b[i]) (complex).
   */
   void expV(DeviceArray<cudaComplex>& a,
             DeviceArray<cudaComplex> const & b,
             const int beginIdA, const int beginIdB, const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _expV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                   b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Elementwise arithmetic square

   /*
   * Vector elementwise square, a[i] = b[i]*b[i] (real).
   */
   void sqV(DeviceArray<cudaReal>& a,
            DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, 
            const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _sqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                  b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector elementwise square, a[i] = b[i]*b[i] (complex).
   */
   void sqV(DeviceArray<cudaComplex>& a,
            DeviceArray<cudaComplex> const & b,
            const int beginIdA, const int beginIdB, 
            const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _sqV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                  b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Elementwise absolute magnitudes

   /*
   * Vector absolute magnitude, a[i] = abs(b[i]) (real).
   */
   void absV(DeviceArray<cudaReal>& a,
             DeviceArray<cudaReal> const & b,
             const int beginIdA, const int beginIdB,
             const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _absV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                   b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector absolute magnitude squared, a[i] = |b[i]|^2 (complex).
   */
   void sqAbsV(DeviceArray<cudaReal>& a,
               DeviceArray<cudaComplex> const & b,
               const int beginIdA, const int beginIdB, 
               const int n)
   {
      UTIL_CHECK(a.capacity() >= n + beginIdA);
      UTIL_CHECK(b.capacity() >= n + beginIdB);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _sqAbsV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA,
                                     b.cArray() + beginIdB, n);
      cudaErrorCheck( cudaGetLastError() );
   }

} // namespace VecOp
} // namespace Pscf
