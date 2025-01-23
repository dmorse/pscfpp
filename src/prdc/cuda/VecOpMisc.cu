/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOpMisc.h"
#include "VecOp.h"
#include <pscf/cuda/ThreadGrid.h>
#include <pscf/cuda/HostDArray.h>
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
   * Vector addition w/ coefficient, a[i] = (b[i]*c) + (d[i]*e), GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array 1 (RHS)
   * \param c  input scalar 1 (RHS)
   * \param d  input array 2 (RHS)
   * \param e  input scalar 2 (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVcVc(cudaReal* a, cudaReal const * b, cudaReal const c, 
                           cudaReal const * d, cudaReal const e, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i] = (b[i] * c) + (d[i] * e);
      }
   }

   /*
   * 3-vector add. w/ coeff, a[i] = (b[i]*c) + (d[i]*e) + (f[i]*g), GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array 1 (RHS)
   * \param c  input scalar 1 (RHS)
   * \param d  input array 2 (RHS)
   * \param e  input scalar 2 (RHS)
   * \param f  input array 3 (RHS)
   * \param g  input scalar 3 (RHS)
   * \param n  size of arrays
   */
   __global__ void _addVcVcVc(cudaReal* a, cudaReal const * b, cudaReal const c, 
                              cudaReal const * d, cudaReal const e, 
                              cudaReal const * f, cudaReal const g, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i] = (b[i] * c) + (d[i] * e) + (f[i] * g);
      }
   }

   /*
   * Vector addition in-place w/ coefficient, a[i] += b[i] * c, GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar
   * \param n  size of arrays
   */
   __global__ void _addEqVc(cudaReal* a, cudaReal const * b, 
                           cudaReal const c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i] += b[i] * c;
      }
   }

   /*
   * Vector subtraction, a[i] = b[i] - c[i] - d, GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input array (RHS)
   * \param d  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _subVVS(cudaReal* a, cudaReal const * b, 
                           cudaReal const * c, cudaReal const d, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i] = b[i] - c[i] - d;
      }
   }

   /*
   * Vector division in-place w/ coeff., a[i] /= (b[i] * c), GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar (RHS)
   * \param n  size of arrays
   */
   __global__ void _divEqVc(cudaComplex* a, cudaReal const * b, 
                           cudaReal const c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i].x /= (b[i] * c);
         a[i].y /= (b[i] * c);
      }
   }

   /*
   * Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param c  input scalar
   * \param n  size of arrays
   */
   __global__ void _expVc(cudaReal* a, cudaReal const * b, 
                        cudaReal const c, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i] = exp(b[i] * c);
      }
   }

   /*
   * Vector assignment in pairs, a1[i] = b[i] and a2[i] = b[i], CUDA kernel.
   *
   * \param a1  output array 1 (LHS)
   * \param a2  output array 2 (LHS)
   * \param s  shared input array to be assigned to both a1 and a2
   * \param n  size of arrays
   */
   __global__ void _eqVPair(cudaReal* a1, cudaReal* a2, 
                           cudaReal const * s, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal input;
      for (int i = startID; i < n; i += nThreads) {
         input = s[i];
         a1[i] = input;
         a2[i] = input;
      }
   }

   /*
   * Vector multiplication in pairs, ax[i] = bx[i] * s[i], CUDA kernel.
   * 
   * \param a1  output array 1 (LHS)
   * \param a2  output array 2 (LHS)
   * \param b1  input array 1 (RHS)
   * \param b2  input array 2 (RHS)
   * \param s  shared input array to be multiplied by both a1 and a2
   * \param n  size of arrays
   */
   __global__ void _mulVVPair(cudaReal* a1, cudaReal * a2, 
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
   * In-place vector multiplication in pairs, ax[i] *= s[i], CUDA kernel
   * 
   * \param a1  output array 1 (LHS)
   * \param a2  output array 2 (LHS)
   * \param s  shared input array to be multiplied by both a1 and a2
   * \param n  size of arrays
   */
   __global__ void _mulEqVPair(cudaReal* a1, cudaReal* a2, 
                              cudaReal const * s, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal input;
      for (int i = startID; i < n; i += nThreads) {
         input = s[i];
         a1[i] *= input;
         a2[i] *= input;
      }
   }

   /*
   * Add an undefined number of vectors pointwise, GPU kernel.
   *
   * The input const pointer 'vecs' points to an array of const pointers. In
   * other words, this is an array of arrays, where each array is represented
   * by its pointer. The size of vecs is nVecs, and the size of vecs[i] is
   * n (if i < nVecs). These nVecs vectors will be added and the result will
   * be stored in vector 'a'.
   *
   * \param a  output array (LHS)
   * \param vecs  array of pointers to DeviceArrays to be added
   * \param nVecs  number of vectors to be added
   * \param n  size of arrays
   */
   __global__ void _addVMany(cudaReal* a, cudaReal const ** vecs,
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
   * Multiply an undefined number of vectors pointwise, GPU kernel.
   *
   * The input const pointer 'vecs' points to an array of const pointers. In
   * other words, this is an array of arrays, where each array is represented
   * by its pointer. The size of vecs is nVecs, and the size of vecs[i] is
   * n (if i < nVecs). These nVecs vectors will be multiplied and the result 
   * will be stored in vector 'a'.
   *
   * \param a  output array (LHS)
   * \param vecs  array of pointers to DeviceArrays to be multiplied
   * \param nVecs  number of vectors to be multiplied
   * \param n  size of arrays
   */
   __global__ void _mulVMany(cudaReal* a, cudaReal const ** vecs,
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
   * Squared norm of complex number, a[i] = norm(b[i])^2, GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _sqNormV(cudaReal* a, cudaComplex const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < n; i += nThreads) {
         a[i] = (b[i].x *  b[i].x) + (b[i].y * b[i].y);
      }
   }

   /*
   * Norm of complex number to the 4th power, a[i] = norm(b[i])^4, GPU kernel.
   *
   * \param a  output array (LHS)
   * \param b  input array (RHS)
   * \param n  size of arrays
   */
   __global__ void _sqSqNormV(cudaReal* a, cudaComplex const * b, const int n)
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal tmp;
      for(int i = startID; i < n; i += nThreads) {
         tmp = (b[i].x *  b[i].x) + (b[i].y * b[i].y);
         a[i] = tmp * tmp;
      }
   }

} // End anonymous namespace


// CUDA kernel wrappers:

// Vector addition w/ coefficient, a[i] = (b[i]*c) + (d[i]*e), kernel wrapper.
void addVcVc(DeviceArray<cudaReal>& a, 
             DeviceArray<cudaReal> const & b, cudaReal const c, 
             DeviceArray<cudaReal> const & d, cudaReal const e)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   UTIL_CHECK(d.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVcVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, d.cArray(), e, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// 3-vector add. w/ coeff, a[i] = (b[i]*c) + (d[i]*e) + (f[i]*g), kernel wrapper
void addVcVcVc(DeviceArray<cudaReal>& a, 
               DeviceArray<cudaReal> const & b, cudaReal const c, 
               DeviceArray<cudaReal> const & d, cudaReal const e, 
               DeviceArray<cudaReal> const & f, cudaReal const g)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   UTIL_CHECK(d.capacity() >= n);
   UTIL_CHECK(f.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVcVcVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, d.cArray(), e, 
                                     f.cArray(), g, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
void addEqVc(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
             cudaReal const c)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector subtraction, a[i] = b[i] - c[i] - d, kernel wrapper.
void subVVS(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
            DeviceArray<cudaReal> const & c, cudaReal const d)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   UTIL_CHECK(c.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), 
                                  c.cArray(), d, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector division in-place w/ coeff., a[i] /= (b[i] * c), kernel wrapper.
void divEqVc(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b, 
             cudaReal const c)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
void expVc(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b, 
           cudaReal const c)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vector assignment in pairs, ax[i] = s[i], kernel wrapper.
void eqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
             DeviceArray<cudaReal> const & s)
{
   const int n = a1.capacity();
   UTIL_CHECK(a2.capacity() == n);
   UTIL_CHECK(s.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqVPair<<<nBlocks, nThreads>>>(a1.cArray(), a2.cArray(), s.cArray(), n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Vec. mul. in pairs, ax[i] = bx[i] * s[i], kernel wrapper.
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
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVVPair<<<nBlocks, nThreads>>>(a1.cArray(), a2.cArray(), b1.cArray(), 
                                     b2.cArray(), s.cArray(), n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// In-place vec. mul. in pairs, ax[i] *= s[i], kernel wrapper
void mulEqVPair(DeviceArray<cudaReal>& a1, DeviceArray<cudaReal>& a2, 
                DeviceArray<cudaReal> const & s)
{
   const int n = a1.capacity();
   UTIL_CHECK(a2.capacity() == n);
   UTIL_CHECK(s.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqVPair<<<nBlocks, nThreads>>>(a1.cArray(), a2.cArray(), s.cArray(), n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Add an undefined number of vectors pointwise, kernel wrapper.
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
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVMany<<<nBlocks, nThreads>>>(a.cArray(), vecs_d.cArray(), nVecs, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Add an undefined number of vectors pointwise, kernel wrapper.
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
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVMany<<<nBlocks, nThreads>>>(a.cArray(), vecs_d.cArray(), nVecs, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Multiply an undefined number of vectors pointwise, kernel wrapper.
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
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVMany<<<nBlocks, nThreads>>>(a.cArray(), vecs_d.cArray(), nVecs, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Multiply an undefined number of vectors pointwise, kernel wrapper.
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
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVMany<<<nBlocks, nThreads>>>(a.cArray(), vecs_d.cArray(), nVecs, n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Squared norm of complex vector, a[i] = norm(b[i])^2, kernel wrapper.
void sqNormV(DeviceArray<cudaReal>& a, DeviceArray<cudaComplex> const & b)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _sqNormV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

// Norm of complex number to the 4th power, a[i] = norm(b[i])^4, kernel wrapper.
void sqSqNormV(DeviceArray<cudaReal>& a, DeviceArray<cudaComplex> const & b)
{
   const int n = a.capacity();
   UTIL_CHECK(b.capacity() >= n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _sqSqNormV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
   cudaErrorCheck( cudaGetLastError() ); // ensure no CUDA errors
}

} // namespace VecOp
} // namespace Cuda
} // namespace Prdc
} // namespace Pscf