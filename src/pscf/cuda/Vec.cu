/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Vec.h"
#include "ThreadGrid.h"
#include <cmath>

namespace Pscf {
namespace Vec {

// Vector assignment, a[i] = b[i], GPU kernel (cudaReal).
__global__ void _eqV(cudaReal* a, cudaReal const * b, const int n) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i];
   }
}

// Vector assignment, a[i] = b[i], GPU kernel (cudaComplex).
__global__ void _eqV(cudaComplex* a, cudaComplex const * b, const int n) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x;
      a[i].y = b[i].y;
   }
}

// Vector assignment, a[i] = b, GPU kernel (cudaReal).
__global__ void _eqS(cudaReal* a, cudaReal const b, const int n) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b;
   }
}

// Vector assignment, a[i] = b, GPU kernel (cudaComplex).
__global__ void _eqS(cudaComplex* a, cudaComplex const b, const int n) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b.x;
      a[i].y = b.y;
   }
}

// Vector assignment, a[i] = b[i], GPU kernel wrapper (cudaReal).
__host__ void eqV(DeviceDArray<cudaReal>& a, 
                  DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector assignment, a[i] = b[i], GPU kernel wrapper (cudaComplex).
__host__ void eqV(DeviceDArray<cudaComplex>& a, 
                  DeviceDArray<cudaComplex> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector assignment, a[i] = b, GPU kernel wrapper (cudaReal).
__host__ void eqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector assignment, a[i] = b, GPU kernel wrapper (cudaComplex).
__host__ void eqS(DeviceDArray<cudaComplex>& a, cudaComplex const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector addition, a[i] = b[i] + c[i], GPU kernel (cudaReal).
__global__ void _addVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] + c[i];
   }
}

// Vector addition, a[i] = b[i] + c[i], GPU kernel (cudaComplex).
__global__ void _addVV(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x + c[i].x;
      a[i].y = b[i].y + c[i].y;
   }
}

// Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed type).
__global__ void _addVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i] + c[i].x;
      a[i].y = c[i].y;
   }
}

// Vector addition, a[i] = b[i] + c, GPU kernel (cudaReal).
__global__ void _addVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] + c;
   }
}

// Vector addition, a[i] = b[i] + c, GPU kernel (cudaComplex).
__global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x + c.x;
      a[i].y = b[i].y + c.y;
   }
}

// Vector addition, a[i] = b[i] + c, GPU kernel (mixed type, b = real).
__global__ void _addVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i] + c.x;
      a[i].y = c.y;
   }
}

// Vector addition, a[i] = b[i] + c, GPU kernel (mixed type, c = real).
__global__ void _addVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x + c;
      a[i].y = b[i].y;
   }
}

// Vector addition, a[i] = b[i] + c[i], GPU kernel wrapper (cudaReal).
__host__ void addVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector addition, a[i] = b[i] + c[i], GPU kernel wrapper (cudaComplex).
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector addition, a[i] = b[i] + c[i], GPU kernel wrapper (mixed type).
__host__ void addVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector addition, a[i] = b[i] + c, GPU kernel wrapper (cudaReal).
__host__ void addVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector addition, a[i] = b[i] + c, GPU kernel wrapper (cudaComplex).
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector addition, a[i] = b[i] + c, GPU kernel wrapper (mixed, b = real).
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    cudaComplex const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector addition, a[i] = b[i] + c, GPU kernel wrapper (mixed, c = real).
__host__ void addVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel (cudaReal).
__global__ void _subVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] - c[i];
   }
}

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel (cudaComplex).
__global__ void _subVV(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x - c[i].x;
      a[i].y = b[i].y - c[i].y;
   }
}

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed type, b = real).
__global__ void _subVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i] - c[i].x;
      a[i].y = 0.0 - c[i].y;
   }
}

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed type, c = real).
__global__ void _subVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x - c[i];
      a[i].y = b[i].y;
   }
}

// Vector subtraction, a[i] = b[i] - c, GPU kernel (cudaReal).
__global__ void _subVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] - c;
   }
}

// Vector subtraction, a[i] = b[i] - c, GPU kernel (cudaComplex).
__global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x - c.x;
      a[i].y = b[i].y - c.y;
   }
}

// Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed type, b = real).
__global__ void _subVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i] - c.x;
      a[i].y = 0.0 - c.y;
   }
}

// Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed type, c = real).
__global__ void _subVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x - c;
      a[i].y = b[i].y;
   }
}

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel wrapper (cudaReal).
__host__ void subVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (mixed type).
__host__ void subVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector subtraction, a[i] = b[i] - c, GPU kernel wrapper (cudaReal).
__host__ void subVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector subtraction, a[i] = b[i] - c, GPU kernel wrapper (cudaComplex).
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    cudaComplex const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
__host__ void subVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector multiplication, a[i] = b[i] * c[i], GPU kernel (cudaReal).
__global__ void _mulVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] * c[i];
   }
}

// Vector multiplication, a[i] = b[i] * c[i], GPU kernel (cudaComplex).
__global__ void _mulVV(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = (b[i].x * c[i].x) - (b[i].y * c[i].y);
      a[i].y = (b[i].x * c[i].y) + (b[i].y * c[i].x);
   }
}

// Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed type).
__global__ void _mulVV(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i] * c[i].x;
      a[i].y = b[i] * c[i].y;
   }
}

// Vector multiplication, a[i] = b[i] * c, GPU kernel (cudaReal).
__global__ void _mulVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] * c;
   }
}

// Vector multiplication, a[i] = b[i] * c, GPU kernel (cudaComplex).
__global__ void _mulVS(cudaComplex* a, cudaComplex const * b, 
                       cudaComplex const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = (b[i].x * c.x) - (b[i].y * c.y);
      a[i].y = (b[i].x * c.y) + (b[i].y * c.x);
   }
}

// Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed type, b = real).
__global__ void _mulVS(cudaComplex* a, cudaReal const * b, 
                       cudaComplex const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i] * c.x;
      a[i].y = b[i] * c.y;
   }
}

// Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed type, c = real).
__global__ void _mulVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x * c;
      a[i].y = b[i].y * c;
   }
}

// Vector multiplication, a[i] = b[i] * c[i], GPU kernel wrapper (cudaReal).
__host__ void mulVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaComplex> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (mixed type).
__host__ void mulVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaComplex> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector multiplication, a[i] = b[i] * c, GPU kernel wrapper (cudaReal).
__host__ void mulVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector multiplication, a[i] = b[i] * c, GPU kernel wrapper (cudaComplex).
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaComplex const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    cudaComplex const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
__host__ void mulVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector division, a[i] = b[i] / c[i], GPU kernel (cudaReal).
__global__ void _divVV(cudaReal* a, cudaReal const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] / c[i];
   }
}

// Vector division, a[i] = b[i] / c[i], GPU kernel (mixed type, c = real).
__global__ void _divVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x / c[i];
      a[i].y = b[i].y / c[i];
   }
}

// Vector division, a[i] = b[i] / c, GPU kernel (cudaReal).
__global__ void _divVS(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = b[i] / c;
   }
}

// Vector division, a[i] = b[i] / c, GPU kernel (mixed type, c = real).
__global__ void _divVS(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x / c;
      a[i].y = b[i].y / c;
   }
}

// Vector division, a[i] = b[i] / c[i], GPU kernel wrapper (cudaReal).
__host__ void divVV(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, 
                    DeviceDArray<cudaReal> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed type, c = real).
__host__ void divVV(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    DeviceDArray<cudaReal> const & c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   UTIL_CHECK(c.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c.cArray(), n);
}

// Vector division, a[i] = b[i] / c, GPU kernel wrapper (cudaReal).
__host__ void divVS(DeviceDArray<cudaReal>& a, 
                    DeviceDArray<cudaReal> const & b, cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
__host__ void divVS(DeviceDArray<cudaComplex>& a, 
                    DeviceDArray<cudaComplex> const & b, 
                    cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Vector exponentiation, a[i] = exp(b[i]), GPU kernel (cudaReal).
__global__ void _expV(cudaReal* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = exp(b[i]);
   }
}

// Vector exponentiation, a[i] = exp(b[i]), GPU kernel (cudaComplex).
__global__ void _expV(cudaComplex* a, cudaComplex const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = exp(b[i].x) * cos(b[i].y);
      a[i].y = exp(b[i].x) * sin(b[i].y);
   }
}

// Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaReal).
__host__ void expV(DeviceDArray<cudaReal>& a, 
                   DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
__host__ void expV(DeviceDArray<cudaComplex>& a, 
                   DeviceDArray<cudaComplex> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector addition in-place, a[i] += b[i], GPU kernel (cudaReal).
__global__ void _addEqV(cudaReal* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] += b[i];
   }
}

// Vector addition in-place, a[i] += b[i], GPU kernel (cudaComplex).
__global__ void _addEqV(cudaComplex* a, cudaComplex const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x += b[i].x;
      a[i].y += b[i].y;
   }
}

// Vector addition in-place, a[i] += b[i], GPU kernel (mixed type).
__global__ void _addEqV(cudaComplex* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x += b[i];
   }
}

// Vector addition in-place, a[i] += b, GPU kernel (cudaReal).
__global__ void _addEqS(cudaReal* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] += b;
   }
}

// Vector addition in-place, a[i] += b, GPU kernel (cudaComplex).
__global__ void _addEqS(cudaComplex* a, cudaComplex const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x += b.x;
      a[i].y += b.y;
   }
}

// Vector addition in-place, a[i] += b, GPU kernel (mixed type).
__global__ void _addEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x += b;
   }
}

// Vector addition in-place, a[i] += b[i], GPU kernel wrapper (cudaReal).
__host__ void
addEqV(DeviceDArray<cudaReal>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector addition in-place, a[i] += b[i], GPU kernel wrapper (cudaComplex).
__host__ void
addEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaComplex> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector addition in-place, a[i] += b[i], GPU kernel wrapper (mixed type).
__host__ void
addEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector addition in-place, a[i] += b, GPU kernel wrapper (cudaReal).
__host__ void
addEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector addition in-place, a[i] += b, GPU kernel wrapper (cudaComplex).
__host__ void
addEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector addition in-place, a[i] += b, GPU kernel wrapper (mixed type).
__host__ void
addEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector subtraction in-place, a[i] -= b[i], GPU kernel (cudaReal).
__global__ void _subEqV(cudaReal* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] -= b[i];
   }
}

// Vector subtraction in-place, a[i] -= b[i], GPU kernel (cudaComplex).
__global__ void _subEqV(cudaComplex* a, cudaComplex const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x -= b[i].x;
      a[i].y -= b[i].y;
   }
}

// Vector subtraction in-place, a[i] -= b[i], GPU kernel (mixed type).
__global__ void _subEqV(cudaComplex* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x -= b[i];
   }
}

// Vector subtraction in-place, a[i] -= b, GPU kernel (cudaReal).
__global__ void _subEqS(cudaReal* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] -= b;
   }
}

// Vector subtraction in-place, a[i] -= b, GPU kernel (cudaComplex).
__global__ void _subEqS(cudaComplex* a, cudaComplex const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x -= b.x;
      a[i].y -= b.y;
   }
}

// Vector subtraction in-place, a[i] -= b, GPU kernel (mixed type).
__global__ void _subEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x -= b;
   }
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
__host__ void
subEqV(DeviceDArray<cudaReal>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
__host__ void
subEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaComplex> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (mixed type).
__host__ void
subEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector subtraction in-place, a[i] -= b, GPU kernel wrapper (cudaReal).
__host__ void
subEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector subtraction in-place, a[i] -= b, GPU kernel wrapper (cudaComplex).
__host__ void
subEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector subtraction in-place, a[i] -= b, GPU kernel wrapper (mixed type).
__host__ void
subEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector multiplication in-place, a[i] *= b[i], GPU kernel (cudaReal).
__global__ void _mulEqV(cudaReal* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] *= b[i];
   }
}

// Vector multiplication in-place, a[i] *= b[i], GPU kernel (cudaComplex).
__global__ void _mulEqV(cudaComplex* a, cudaComplex const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   cudaComplex c;
   for(int i = startID; i < n; i += nThreads) {
      c.x = (a[i].x * b[i].x) - (a[i].y * b[i].y);
      c.y = (a[i].x * b[i].y) + (a[i].y * b[i].x);
      a[i].x = c.x;
      a[i].y = c.y;
   }
}

// Vector multiplication in-place, a[i] *= b[i], GPU kernel (mixed type).
__global__ void _mulEqV(cudaComplex* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x *= b[i];
      a[i].y *= b[i];
   }
}

// Vector multiplication in-place, a[i] *= b, GPU kernel (cudaReal).
__global__ void _mulEqS(cudaReal* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] *= b;
   }
}

// Vector multiplication in-place, a[i] *= b, GPU kernel (cudaComplex).
__global__ void _mulEqS(cudaComplex* a, cudaComplex const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   cudaComplex c;
   for(int i = startID; i < n; i += nThreads) {
      c.x = (a[i].x * b.x) - (a[i].y * b.y);
      c.y = (a[i].x * b.y) + (a[i].y * b.x);
      a[i].x = c.x;
      a[i].y = c.y;
   }
}

// Vector multiplication in-place, a[i] *= b, GPU kernel (mixed type).
__global__ void _mulEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x *= b;
      a[i].y *= b;
   }
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaReal).
__host__ void
mulEqV(DeviceDArray<cudaReal>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaComplex).
__host__ void
mulEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaComplex> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (mixed type).
__host__ void
mulEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
__host__ void
mulEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
__host__ void
mulEqS(DeviceDArray<cudaComplex>& a, cudaComplex const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (mixed type).
__host__ void
mulEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector division in-place, a[i] /= b[i], GPU kernel (cudaReal).
__global__ void _divEqV(cudaReal* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] /= b[i];
   }
}

// Vector division in-place, a[i] /= b[i], GPU kernel (mixed type).
__global__ void _divEqV(cudaComplex* a, cudaReal const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x /= b[i];
      a[i].y /= b[i];
   }
}

// Vector division in-place, a[i] /= b, GPU kernel (cudaReal).
__global__ void _divEqS(cudaReal* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] /= b;
   }
}

// Vector division in-place, a[i] /= b, GPU kernel (mixed type).
__global__ void _divEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x /= b;
      a[i].y /= b;
   }
}

// Vector division in-place, a[i] /= b[i], GPU kernel wrapper (cudaReal).
__host__ void
divEqV(DeviceDArray<cudaReal>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector division in-place, a[i] /= b[i], GPU kernel wrapper (mixed type).
__host__ void
divEqV(DeviceDArray<cudaComplex>& a, DeviceDArray<cudaReal> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector division in-place, a[i] /= b, GPU kernel wrapper (cudaReal).
__host__ void
divEqS(DeviceDArray<cudaReal>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector division in-place, a[i] /= b, GPU kernel wrapper (mixed type).
__host__ void
divEqS(DeviceDArray<cudaComplex>& a, cudaReal const b)
{
   int n = a.capacity();
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqS<<<nBlocks, nThreads>>>(a.cArray(), b, n);
}

// Vector addition in-place w/ coefficient, a[i] += b[i] * c, GPU kernel.
__global__ void _addEqMulVS(cudaReal* a, cudaReal const * b, 
                            cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] += b[i] * c;
   }
}

// Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
__host__ void addEqMulVS(DeviceDArray<cudaReal>& a, 
                         DeviceDArray<cudaReal> const & b, 
                         cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqMulVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

// Squared norm of complex vector, a[i] = norm(b[i])^2, GPU kernel.
__global__ void _sqNormV(cudaReal* a, cudaComplex const * b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = (b[i].x *  b[i].x) + (b[i].y * b[i].y);
   }
}

// Squared norm of complex vector, a[i] = norm(b[i])^2, kernel wrapper.
__host__ void sqNormV(DeviceDArray<cudaReal>& a, 
                      DeviceDArray<cudaComplex> const & b)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _sqNormV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), n);
}

// Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), GPU kernel.
__global__ void _expMulVS(cudaReal* a, cudaReal const * b, 
                          cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = exp(b[i] * c);
   }
}

// Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
__host__ void expMulVS(DeviceDArray<cudaReal>& a, 
                       DeviceDArray<cudaReal> const & b, 
                       cudaReal const c)
{
   int n = a.capacity();
   UTIL_CHECK(b.capacity() == n);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expMulVS<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), c, n);
}

}
}
