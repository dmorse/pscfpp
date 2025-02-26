/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOp.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/cudaErrorCheck.h>
#include <cmath>

namespace Pscf {
namespace Prdc {
namespace Cuda {
namespace VecOp {

// CUDA kernels:
// (defined in anonymous namespace, used within this file only)

namespace {

   /*
   * Vector assignment, a[i] = b[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _eqV(cudaReal* a, cudaReal const * b, const int n) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b[i];
      }
   }

   /*
   * Vector assignment, a[i] = b[i], GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _eqV(cudaComplex* a, cudaComplex const * b, const int n) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x = b[i].x;
         a[i].y = b[i].y;
      }
   }

   /*
   * Vector assignment, a[i] = b, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _eqS(cudaReal* a, const cudaReal b, const int n) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b;
      }
   }

   /*
   * Vector assignment, a[i] = b, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _eqS(cudaComplex* a, const cudaComplex b, const int n) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x = b.x;
         a[i].y = b.y;
      }
   }

   /*
   * Vector addition, a[i] = b[i] + c[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVV(cudaReal* a, cudaReal const * b, 
                          cudaReal const * c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b[i] + c[i];
      }
   }

   /*
   * Vector addition, a[i] = b[i] + c[i], GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVV(cudaComplex* a, cudaComplex const * b, 
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
   * Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVV(cudaComplex* a, cudaReal const * b, 
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
   * Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, c = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVV(cudaComplex* a, cudaComplex const * b, 
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
   * Vector addition, a[i] = b[i] + c, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVS(cudaReal* a, cudaReal const * b, 
                          const cudaReal c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b[i] + c;
      }
   }

   /*
   * Vector addition, a[i] = b[i] + c, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
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
   * Vector addition, a[i] = b[i] + c, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVS(cudaComplex* a, cudaReal const * b, 
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
   * Vector addition, a[i] = b[i] + c, GPU kernel (mixed, c = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
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
   * Vector subtraction, a[i] = b[i] - c[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVV(cudaReal* a, cudaReal const * b, 
                          cudaReal const * c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b[i] - c[i];
      }
   }

   /*
   * Vector subtraction, a[i] = b[i] - c[i], GPU kernel (cudaComplex).
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
   * Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, b = real).
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
   * Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, c = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVV(cudaComplex* a, cudaComplex const * b, 
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
   * Vector subtraction, a[i] = b[i] - c, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVS(cudaReal* a, cudaReal const * b, 
                          const cudaReal c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b[i] - c;
      }
   }

   /*
   * Vector subtraction, a[i] = b[i] - c, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
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
   * Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVS(cudaComplex* a, cudaReal const * b, 
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
   * Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, c = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
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
   * Vector multiplication, a[i] = b[i] * c[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulVV(cudaReal* a, cudaReal const * b, 
                          cudaReal const * c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b[i] * c[i];
      }
   }

   /*
   * Vector multiplication, a[i] = b[i] * c[i], GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulVV(cudaComplex* a, cudaComplex const * b, 
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
   * Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulVV(cudaComplex* a, cudaReal const * b, 
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
   * Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, c = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulVV(cudaComplex* a, cudaComplex const * b, 
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
   * Vector multiplication, a[i] = b[i] * c, GPU kernel (cudaReal).
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
   * Vector multiplication, a[i] = b[i] * c, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulVS(cudaComplex* a, cudaComplex const * b, 
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
   * Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulVS(cudaComplex* a, cudaReal const * b, 
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
   * Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, c = real).
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
   * Vector division, a[i] = b[i] / c[i], GPU kernel (cudaReal).
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
   * Vector division, a[i] = b[i] / c[i], GPU kernel (mixed, c = real).
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
   * Vector division, a[i] = b[i] / c, GPU kernel (cudaReal).
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
   * Vector division, a[i] = b[i] / c, GPU kernel (mixed, c = real).
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
   * Vector division, a[i] = b / c[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param c  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _divSV(cudaReal* a, 
                          const cudaReal b, 
                          cudaReal const * c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = b / c[i];
      }
   }

   /*
   * Vector exponentiation, a[i] = exp(b[i]), GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _expV(cudaReal* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] = exp(b[i]);
      }
   }

   /*
   * Vector exponentiation, a[i] = exp(b[i]), GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _expV(cudaComplex* a, cudaComplex const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x = exp(b[i].x) * cos(b[i].y);
         a[i].y = exp(b[i].x) * sin(b[i].y);
      }
   }

   /*
   * Vector addition in-place, a[i] += b[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addEqV(cudaReal* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] += b[i];
      }
   }

   /*
   * Vector addition in-place, a[i] += b[i], GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addEqV(cudaComplex* a, cudaComplex const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x += b[i].x;
         a[i].y += b[i].y;
      }
   }

   /*
   * Vector addition in-place, a[i] += b[i], GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _addEqV(cudaComplex* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x += b[i];
      }
   }

   /*
   * Vector addition in-place, a[i] += b, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addEqS(cudaReal* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] += b;
      }
   }

   /*
   * Vector addition in-place, a[i] += b, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addEqS(cudaComplex* a, const cudaComplex b, 
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
   * Vector addition in-place, a[i] += b, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _addEqS(cudaComplex* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x += b;
      }
   }

   /*
   * Vector subtraction in-place, a[i] -= b[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _subEqV(cudaReal* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] -= b[i];
      }
   }

   /*
   * Vector subtraction in-place, a[i] -= b[i], GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _subEqV(cudaComplex* a, cudaComplex const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x -= b[i].x;
         a[i].y -= b[i].y;
      }
   }

   /*
   * Vector subtraction in-place, a[i] -= b[i], GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _subEqV(cudaComplex* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x -= b[i];
      }
   }

   /*
   * Vector subtraction in-place, a[i] -= b, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subEqS(cudaReal* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] -= b;
      }
   }

   /*
   * Vector subtraction in-place, a[i] -= b, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subEqS(cudaComplex* a, const cudaComplex b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x -= b.x;
         a[i].y -= b.y;
      }
   }

   /*
   * Vector subtraction in-place, a[i] -= b, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subEqS(cudaComplex* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x -= b;
      }
   }

   /*
   * Vector multiplication in-place, a[i] *= b[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulEqV(cudaReal* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] *= b[i];
      }
   }

   /*
   * Vector multiplication in-place, a[i] *= b[i], GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulEqV(cudaComplex* a, cudaComplex const * b, const int n)
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
   * Vector multiplication in-place, a[i]*=b[i], GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulEqV(cudaComplex* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x *= b[i];
         a[i].y *= b[i];
      }
   }

   /*
   * Vector multiplication in-place, a[i] *= b, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulEqS(cudaReal* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] *= b;
      }
   }

   /*
   * Vector multiplication in-place, a[i] *= b, GPU kernel (cudaComplex).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulEqS(cudaComplex* a, const cudaComplex b, const int n)
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
   * Vector multiplication in-place, a[i] *= b, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _mulEqS(cudaComplex* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x *= b;
         a[i].y *= b;
      }
   }

   /*
   * Vector division in-place, a[i] /= b[i], GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _divEqV(cudaReal* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] /= b[i];
      }
   }

   /*
   * Vector division in-place, a[i] /= b[i], GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _divEqV(cudaComplex* a, cudaReal const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x /= b[i];
         a[i].y /= b[i];
      }
   }

   /*
   * Vector division in-place, a[i] /= b, GPU kernel (cudaReal).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _divEqS(cudaReal* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i] /= b;
      }
   }

   /*
   * Vector division in-place, a[i] /= b, GPU kernel (mixed, b = real).
   *
   * \param a  output array (LHS)
   * \param b  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _divEqS(cudaComplex* a, const cudaReal b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < n; i += nThreads) {
         a[i].x /= b;
         a[i].y /= b;
      }
   }

} // end anonymous namespace


// CUDA kernel wrappers:

// Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
void eqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector assignment, a[i] = b[i], kernel wrapper (cudaComplex).
void eqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector assignment, a[i] = b, kernel wrapper (cudaReal).
void eqS(DeviceArray<cudaReal>& a, const cudaReal b,
         const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
void eqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
         const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaReal).
void addVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
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
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaComplex).
void addVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, b = real).
void addVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, c = real).
void addVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
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
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (cudaReal).
void addVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (cudaComplex).
void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           const cudaComplex c, const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, b = real).
void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           const cudaComplex c, const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, c = real).
void addVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaReal).
void subVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
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
   _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n)
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i]=b[i]-c[i], kernel wrapper (mixed, b=real).
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n)
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i]=b[i]-c[i], kernel wrapper (mixed, c=real).
void subVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
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
   _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaReal).
void subVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, const cudaReal c,
                    const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaComplex).
void subVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
void subVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
void subVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaReal).
void mulVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
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
   _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n)
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, b = real).
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
           DeviceArray<cudaComplex> const & c, const int beginIdA, 
           const int beginIdB, const int beginIdC, const int n)
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, c = real).
void mulVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
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
   _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
void mulVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
void mulVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           const cudaComplex c, const int beginIdA, const int beginIdB, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
void mulVS(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division, a[i] = b[i] / c[i], kernel wrapper (cudaReal).
void divVV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed, c = real).
void divVV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b, 
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division, a[i] = b[i] / c, kernel wrapper (cudaReal).
void divVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB, 
           const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
void divVS(DeviceArray<cudaComplex>& a, 
           DeviceArray<cudaComplex> const & b, 
           const cudaReal c, const int beginIdA, const int beginIdB, 
           const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Division of scalar by  vector, a[i] = b / c[i], kernel wrapper (cudaReal).
void divSV(DeviceArray<cudaReal>& a, 
           const cudaReal b, 
           DeviceArray<cudaReal> const & c, 
           const int beginIdA, const int beginIdC, 
           const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divSV<<<nBlocks, nThreads>>>(a.cArray() + beginIdA, b, 
                                 c.cArray() + beginIdC, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
void expV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
          const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                b.cArray()+beginIdB, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
void expV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
          const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                b.cArray()+beginIdB, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
void addEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place, a[i] += b[i], kernel wrapper (cudaComplex).
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place, a[i] += b[i], kernel wrapper (mixed).
void addEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
void addEqS(DeviceArray<cudaReal>& a, const cudaReal b, 
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
void addEqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place, a[i] += b, kernel wrapper (mixed).
void addEqS(DeviceArray<cudaComplex>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (mixed).
void subEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
void subEqS(DeviceArray<cudaReal>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
void subEqS(DeviceArray<cudaComplex>& a, const cudaComplex b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed).
void subEqS(DeviceArray<cudaComplex>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
void mulEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaComplex).
void mulEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaComplex> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (mixed).
void mulEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
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
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
void mulEqS(DeviceArray<cudaReal>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
void mulEqS(DeviceArray<cudaComplex>& a, 
            const cudaComplex b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (mixed).
void mulEqS(DeviceArray<cudaComplex>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
void divEqV(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division in-place, a[i] /= b[i], kernel wrapper (mixed).
void divEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
            const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
void divEqS(DeviceArray<cudaReal>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division in-place, a[i] /= b, kernel wrapper (mixed).
void divEqS(DeviceArray<cudaComplex>& a, const cudaReal b,
            const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

} // namespace VecOp
} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
