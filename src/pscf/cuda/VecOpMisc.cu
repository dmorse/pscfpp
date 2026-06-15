/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOpMisc.h"
#include "VecOp.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaErrorCheck.h>
#include <cmath>

namespace Pscf {
namespace VecOp {

   // Anonymous namespace for CUDA kernels (only accessible in this file)
   namespace {

      /*
      * Add scaled vectors, a[i] = (b[i]*c) + (d[i]*e).
      *
      * \param a  real output array (LHS)
      * \param b1  real input array 1 (RHS)
      * \param c1  real coefficient of array 1 (RHS)
      * \param b2  real input array 2 (RHS)
      * \param c2  real coefficient of array 2 (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVcVc(cudaReal* a,
                    cudaReal const * b1, cudaReal const c1,
                    cudaReal const * b2, cudaReal const c2,
                    const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = (b1[i] * c1) + (b2[i] * c2);
         }
      }

      /*
      * Add a scaled vector and a scalar, a[i] = b[i]*c + s.
      *
      * \param a  output array (LHS)
      * \param b  real input array (RHS)
      * \param c  real coefficient of array b (RHS)
      * \param s  real additive scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVcS(cudaReal* a,
                   cudaReal const * b, cudaReal const c,
                   cudaReal const s,
                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] * c + s;
         }
      }

      /*
      * Add scaled vector in-place, a[i] += b[i] * c.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar
      * \param n  size of arrays
      */
      __global__
      void _addEqVc(cudaReal* a,
                    cudaReal const * b, cudaReal const c,
                    const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b[i] * c;
         }
      }

      /*
      * Add 3 scaled vectors, a[i] = b[i]*c + d[i]*e + f[i]*g.
      *
      * \param a  output array (LHS)
      * \param b1  input array 1 (RHS)
      * \param c1  input scalar 1 (RHS)
      * \param b2  input array 2 (RHS)
      * \param c2  input scalar 2 (RHS)
      * \param b3  input array 3 (RHS)
      * \param c3  input scalar 3 (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVcVcVc(cudaReal* a,
                      cudaReal const * b1, cudaReal const c1,
                      cudaReal const * b2, cudaReal const c2,
                      cudaReal const * b3, cudaReal const c3,
                      const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = (b1[i] * c1) + (b2[i] * c2) + (b3[i] * c3);
         }
      }

      /*
      * Add scaled vectors and scalar, a[i] = b1[i]*c1 + b2[i]*c2 + s.
      *
      * \param a  real output array (LHS)
      * \param b1  real input array 1 (RHS)
      * \param c1  real coefficient of b1 (RHS)
      * \param b2  real input array 2 (RHS)
      * \param c2  real coefficient of b2 (RHS)
      * \param s  real scalar summand (RHS)
      * \param n  size of arrays
      */
      __global__
      void _addVcVcS(cudaReal* a,
                    cudaReal const * b1, cudaReal const c1,
                    cudaReal const * b2, cudaReal const c2,
                    const cudaReal s, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = (b1[i] * c1) + (b2[i] * c2) + s;
         }
      }

      // In-place addition of scaled vectors

      #if 0
      /*
      * Vector subtraction, a[i] = b[i] - c[i] - d.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input array (RHS)
      * \param d  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _subVVS(cudaReal* a,
                   cudaReal const * b,
                   cudaReal const * c,
                   cudaReal const d,
                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = b[i] - c[i] - d;
         }
      }
      #endif

      /*
      * Vector division in-place w/ coeff., a[i] /= (b[i] * c).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar (RHS)
      * \param n  size of arrays
      */
      __global__
      void _divEqVc(cudaComplex* a,
                    cudaReal const * b,
                    cudaReal const c,
                    const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i].x /= (b[i] * c);
            a[i].y /= (b[i] * c);
         }
      }

      /*
      * Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c).
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param c  input scalar
      * \param n  size of arrays
      */
      __global__
      void _expVc(cudaReal* a,
                  cudaReal const * b,
                  cudaReal const c,
                  const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = exp(b[i] * c);
         }
      }

      // Pair operations

      /*
      * Vector assignment in pairs, a1[i] = b[i] and a2[i] = b[i].
      *
      * \param a1  output array 1 (LHS)
      * \param a2  output array 2 (LHS)
      * \param b  shared input array assigned to both a1 and a2 (RHS)
      * \param n  size of arrays
      */
      __global__
      void _eqVPair(cudaReal* a1, cudaReal* a2,
                    cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal input;
         for (int i = startID; i < n; i += nThreads) {
            input = b[i];
            a1[i] = input;
            a2[i] = input;
         }
      }

      /*
      * Vector multiplication in pairs, ax[i] = bx[i] * s[i].
      *
      * \param a1  output array 1 (LHS)
      * \param a2  output array 2 (LHS)
      * \param b1  input array 1 (RHS)
      * \param b2  input array 2 (RHS)
      * \param s  shared input array to be multiplied by both a1 and a2
      * \param n  size of arrays
      */
      __global__
      void _mulVVPair(cudaReal* a1, cudaReal * a2,
                      cudaReal const * b1, cudaReal const * b2,
                      cudaReal const * s, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal input;
         for (int i = startID; i < n; i += nThreads) {
            input = s[i];
            a1[i] = b1[i] * input;
            a2[i] = b2[i] * input;
         }
      }

      /*
      * In-place vector-scalar multiplication in pairs, ax[i] *= s[i].
      *
      * \param a1  output array 1 (LHS)
      * \param a2  output array 2 (LHS)
      * \param b  shared input array to multiplied both a1 and a2
      * \param n  size of arrays
      */
      __global__
      void _mulEqVPair(cudaReal* a1, cudaReal* a2,
                       cudaReal const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal input;
         for (int i = startID; i < n; i += nThreads) {
            input = b[i];
            a1[i] *= input;
            a2[i] *= input;
         }
      }

      /*
      * Add an undefined number of vectors pointwise.
      *
      * The input const pointer 'vecs' points to an array of const pointers.
      * In other words, this is an array of arrays, where each array is
      * represented by its pointer. The size of vecs is nVecs, and the size
      * of vecs[i] is n (if i < nVecs). These nVecs vectors will be added
      * and the result will be stored in vector 'a'.
      *
      * \param a  output array (LHS)
      * \param vecs  array of pointers to DeviceArrays to be added
      * \param nVecs  number of vectors to be added
      * \param n  size of arrays
      */
      __global__
      void _addVMany(cudaReal* a,
                     cudaReal const ** vecs,
                     const int nVecs, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            cudaReal sum = vecs[0][i];
            for (int j = 1; j < nVecs; j++) {
               sum += vecs[j][i];
            }
            a[i] = sum;
         }
      }

      /*
      * Multiply an undefined number of vectors pointwise.
      *
      * The input pointer 'vecs' points to an array of const pointers.
      * In other words, this is an array of arrays, where each array is
      * represented by its pointer. The size of vecs is nVecs, and the
      * size of vecs[i] is n (if i < nVecs). These nVecs vectors will
      * be multiplied and the result will be stored in vector 'a'.
      *
      * \param a  output array (LHS)
      * \param vecs  array of pointers to DeviceArrays to be multiplied
      * \param nVecs  number of vectors to be multiplied
      * \param n  size of arrays
      */
      __global__
      void _mulVMany(cudaReal* a,
                     cudaReal const ** vecs,
                     const int nVecs, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            cudaReal prod = vecs[0][i];
            for (int j = 1; j < nVecs; j++) {
               prod *= vecs[j][i];
            }
            a[i] = prod;
         }
      }

      /*
      * Fourth power of magnitude of complex number, a[i] = |b[i]|^4.
      *
      * \param a  output array (LHS)
      * \param b  input array (RHS)
      * \param n  size of arrays
      */
      __global__
      void _sqSqAbsV(cudaReal* a, cudaComplex const * b, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal tmp;
         for (int i = startID; i < n; i += nThreads) {
            tmp = (b[i].x *  b[i].x) + (b[i].y * b[i].y);
            a[i] = tmp * tmp;
         }
      }

   } // End anonymous namespace

   // Linear combinations - addition of scaled vectors

   /*
   * Add two scaled vectors, a[i] = b1[i]*c1 + b2[i]*c2.
   */
   void addVcVc(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b1, cudaReal const c1,
                DeviceArray<cudaReal> const & b2, cudaReal const c2)
   {
      const int n = a.capacity();
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVcVc<<<nBlocks, nThreads>>>(a.cArray(),
                                      b1.cArray(), c1,
                                      b2.cArray(), c2, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Add a scaled vector and a scalar, a[i] = b[i]*c + s.
   */
   void addVcS(DeviceArray<cudaReal>& a,
               DeviceArray<cudaReal> const & b, cudaReal const c,
               cudaReal const s)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVcS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, s, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Add 3 scaled vectors, a[i] = b1[i]*c1 + b2[i]*c2 + b3[i]*c3).
   */
   void addVcVcVc(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b1, cudaReal const c1,
                  DeviceArray<cudaReal> const & b2, cudaReal const c2,
                  DeviceArray<cudaReal> const & b3, cudaReal const c3)
   {
      const int n = a.capacity();
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);
      UTIL_CHECK(b3.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVcVcVc<<<nBlocks, nThreads>>>(a.cArray(),
                                        b1.cArray(), c1,
                                        b2.cArray(), c2,
                                        b3.cArray(), c3, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Add 2 scaled vectors and scalar, a[i] = b1[i]*c1 + b2[i]*c2 + s;
   */
   void addVcVcS(DeviceArray<cudaReal>& a,
                  DeviceArray<cudaReal> const & b1, cudaReal const c1,
                  DeviceArray<cudaReal> const & b2, cudaReal const c2,
                  cudaReal const s)
   {
      const int n = a.capacity();
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVcVcS<<<nBlocks, nThreads>>>(a.cArray(),
                                       b1.cArray(), c1,
                                       b2.cArray(), c2,
                                       s, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   // Compound (in-place) addition of scaled vectors.

   /*
   * Add scaled vector in-place, a[i] += b[i] * c.
   */
   void addEqVc(DeviceArray<cudaReal>& a,
                DeviceArray<cudaReal> const & b,
                cudaReal const c)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector in-place division by scaled vector, a[i] /= (b[i] * c).
   */
   void divEqVc(DeviceArray<cudaComplex>& a,
                DeviceArray<cudaReal> const & b, cudaReal const c)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _divEqVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c).
   */
   void expVc(DeviceArray<cudaReal>& a,
              DeviceArray<cudaReal> const & b, cudaReal const c)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _expVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vector assignment in pairs, ax[i] = s[i].
   */
   void eqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2,
                DeviceArray<cudaReal> const & s)
   {
      const int n = a1.capacity();
      UTIL_CHECK(a2.capacity() == n);
      UTIL_CHECK(s.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _eqVPair<<<nBlocks, nThreads>>>(a1.cArray(), a2.cArray(),
                                      s.cArray(), n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Vec. mul. in pairs, ax[i] = bx[i] * s[i].
   */
   void mulVVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2,
                  DeviceArray<cudaReal> const & b1,
                  DeviceArray<cudaReal> const & b2,
                  DeviceArray<cudaReal> const & s)
   {
      const int n = a1.capacity();
      UTIL_CHECK(a2.capacity() == n);
      UTIL_CHECK(b1.capacity() >= n);
      UTIL_CHECK(b2.capacity() >= n);
      UTIL_CHECK(s.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVVPair<<<nBlocks, nThreads>>>(a1.cArray(), a2.cArray(),
                                        b1.cArray(), b2.cArray(),
                                        s.cArray(), n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * In-place vec. mul. in pairs, ax[i] *= s[i].
   */
   void mulEqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2,
                   DeviceArray<cudaReal> const & s)
   {
      const int n = a1.capacity();
      UTIL_CHECK(a2.capacity() == n);
      UTIL_CHECK(s.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulEqVPair<<<nBlocks, nThreads>>>(a1.cArray(), a2.cArray(),
                                         s.cArray(), n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Add an undefined number of vectors pointwise.
   */
   void addVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> > const & vecs)
   {
      int nVecs = vecs.capacity();
      UTIL_CHECK(nVecs > 1);
      int n = vecs[0].capacity();

      if (nVecs == 2) {
         addVV(a, vecs[0], vecs[1]);
         return;
      }

      // Create array of pointers to arrays on host
      HostDArray<cudaReal const *> vecs_h(nVecs);
      for (int i = 0; i < nVecs; i++) {
         vecs_h[i] = vecs[i].cArray();
      }
      DeviceArray<cudaReal const *> vecs_d(nVecs);
      vecs_d = vecs_h; // transfer array of pointers to device

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVMany<<<nBlocks, nThreads>>>(a.cArray(),
                                       vecs_d.cArray(), nVecs, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Add an undefined number of vectors pointwise.
   */
   void addVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> const *> const & vecs)
   {
      int nVecs = vecs.capacity();
      UTIL_CHECK(nVecs > 1);
      int n = vecs[0]->capacity();

      if (nVecs == 2) {
         addVV(a, *vecs[0], *vecs[1]);
         return;
      }

      // Create array of pointers to arrays on host
      HostDArray<cudaReal const *> vecs_h(nVecs);
      for (int i = 0; i < nVecs; i++) {
         vecs_h[i] = vecs[i]->cArray();
      }
      DeviceArray<cudaReal const *> vecs_d(nVecs);
      vecs_d = vecs_h; // transfer array of pointers to device

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addVMany<<<nBlocks, nThreads>>>(a.cArray(),
                                       vecs_d.cArray(), nVecs, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Multiply an undefined number of vectors pointwise.
   */
   void mulVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> > const & vecs)
   {
      int nVecs = vecs.capacity();
      UTIL_CHECK(nVecs > 1);
      int n = vecs[0].capacity();

      if (nVecs == 2) {
         mulVV(a, vecs[0], vecs[1]);
         return;
      }

      // Create array of pointers to arrays on host
      HostDArray<cudaReal const *> vecs_h(nVecs);
      for (int i = 0; i < nVecs; i++) {
         vecs_h[i] = vecs[i].cArray();
      }
      DeviceArray<cudaReal const *> vecs_d(nVecs);
      vecs_d = vecs_h; // transfer array of pointers to device

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVMany<<<nBlocks, nThreads>>>(a.cArray(), vecs_d.cArray(),
		                       nVecs, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Multiply an undefined number of vectors pointwise.
   */
   void mulVMany(DeviceArray<cudaReal>& a,
                 DArray<DeviceArray<cudaReal> const *> const & vecs)
   {
      int nVecs = vecs.capacity();
      UTIL_CHECK(nVecs > 1);
      int n = vecs[0]->capacity();

      if (nVecs == 2) {
         mulVV(a, *vecs[0], *vecs[1]);
         return;
      }

      // Create array of pointers to arrays on host
      HostDArray<cudaReal const *> vecs_h(nVecs);
      for (int i = 0; i < nVecs; i++) {
         vecs_h[i] = vecs[i]->cArray();
      }
      DeviceArray<cudaReal const *> vecs_d(nVecs);
      vecs_d = vecs_h; // transfer array of pointers to device

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _mulVMany<<<nBlocks, nThreads>>>(a.cArray(),
                                       vecs_d.cArray(), nVecs, n);
      cudaErrorCheck( cudaGetLastError() );
   }

   /*
   * Fourth power of magnitude of complex vector, a[i] = |b[i]|^4.
   */
   void sqSqAbsV(DeviceArray<cudaReal>& a,
                DeviceArray<cudaComplex> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _sqSqAbsV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
      cudaErrorCheck( cudaGetLastError() );
   }

} // namespace VecOp
} // namespace Pscf
