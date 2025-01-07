/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOp.h"
#include "ThreadGrid.h"
#include "HostDArray.tpp"  // tpp needed to use implicit instantiation
#include "DeviceArray.tpp" // tpp needed to use implicit instantiation
#include <cmath>

namespace Pscf {
namespace VecOp {

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

// Vector assignment, a[i] = b[i], kernel wrapper (cudaReal).
__host__ void eqV(DeviceArray<cudaReal>& a, 
                  DeviceArray<cudaReal> const & b,
                  const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                               b.cArray()+beginIdB, n);
}

// Vector assignment, a[i] = b[i], kernel wrapper (cudaComplex).
__host__ void eqV(DeviceArray<cudaComplex>& a, 
                  DeviceArray<cudaComplex> const & b,
                  const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                               b.cArray()+beginIdB, n);
}

// Vector assignment, a[i] = b, kernel wrapper (cudaReal).
__host__ void eqS(DeviceArray<cudaReal>& a, cudaReal const b,
                  const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector assignment, a[i] = b, kernel wrapper (cudaComplex).
__host__ void eqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                  const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _eqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
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

// Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, b = real).
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

// Vector addition, a[i] = b[i] + c[i], GPU kernel (mixed, c = real).
__global__ void _addVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x + c[i];
      a[i].y = b[i].y;
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

// Vector addition, a[i] = b[i] + c, GPU kernel (mixed, b = real).
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

// Vector addition, a[i] = b[i] + c, GPU kernel (mixed, c = real).
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

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaReal).
__host__ void addVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (cudaComplex).
__host__ void addVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, b = real).
__host__ void addVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector addition, a[i] = b[i] + c[i], kernel wrapper (mixed, c = real).
__host__ void addVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (cudaReal).
__host__ void addVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (cudaComplex).
__host__ void addVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaComplex const c,
                    const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, b = real).
__host__ void addVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    cudaComplex const c,
                    const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector addition, a[i] = b[i] + c, kernel wrapper (mixed, c = real).
__host__ void addVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaReal const c,
                    const int beginIdA, const int beginIdB,int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
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

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, b = real).
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

// Vector subtraction, a[i] = b[i] - c[i], GPU kernel (mixed, c = real).
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

// Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, b = real).
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

// Vector subtraction, a[i] = b[i] - c, GPU kernel (mixed, c = real).
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

// Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaReal).
__host__ void subVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector subtraction, a[i] = b[i] - c[i], kernel wrapper (cudaComplex).
__host__ void subVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector subtraction, a[i]=b[i]-c[i], kernel wrapper (mixed, b=real).
__host__ void subVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector subtraction, a[i]=b[i]-c[i], kernel wrapper (mixed, c=real).
__host__ void subVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaReal).
__host__ void subVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (cudaComplex).
__host__ void subVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, b = real).
__host__ void subVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector subtraction, a[i] = b[i] - c, kernel wrapper (mixed, c = real).
__host__ void subVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaReal const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
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

// Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, b = real).
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

// Vector multiplication, a[i] = b[i] * c[i], GPU kernel (mixed, c = real).
__global__ void _mulVV(cudaComplex* a, cudaComplex const * b, 
                       cudaReal const * c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x = b[i].x * c[i];
      a[i].y = b[i].y * c[i];
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

// Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, b = real).
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

// Vector multiplication, a[i] = b[i] * c, GPU kernel (mixed, c = real).
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

// Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaReal).
__host__ void mulVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector multiplication, a[i] = b[i] * c[i], kernel wrapper (cudaComplex).
__host__ void mulVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, b = real).
__host__ void mulVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaComplex> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector multiplication, a[i]=b[i]*c[i], kernel wrapper (mixed, c = real).
__host__ void mulVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaReal).
__host__ void mulVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (cudaComplex).
__host__ void mulVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, b = real).
__host__ void mulVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaReal> const & b, 
                    cudaComplex const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector multiplication, a[i] = b[i] * c, kernel wrapper (mixed, c = real).
__host__ void mulVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaReal const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
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

// Vector division, a[i] = b[i] / c[i], GPU kernel (mixed, c = real).
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

// Vector division, a[i] = b[i] / c, GPU kernel (mixed, c = real).
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

// Vector division, a[i] = b[i] / c[i], kernel wrapper (cudaReal).
__host__ void divVV(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector division, a[i] = b[i] / c[i], kernel wrapper (mixed, c = real).
__host__ void divVV(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    DeviceArray<cudaReal> const & c,
                    const int beginIdA, const int beginIdB, 
                    const int beginIdC, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   UTIL_CHECK(c.capacity() >= n + beginIdC);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c.cArray()+beginIdC, n);
}

// Vector division, a[i] = b[i] / c, kernel wrapper (cudaReal).
__host__ void divVS(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, cudaReal const c,
                    const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
}

// Vector division, a[i] = b[i] / c, kernel wrapper (mixed, c = real).
__host__ void divVS(DeviceArray<cudaComplex>& a, 
                    DeviceArray<cudaComplex> const & b, 
                    cudaReal const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divVS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b.cArray()+beginIdB, 
                                 c, n);
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
__host__ void expV(DeviceArray<cudaReal>& a, 
                   DeviceArray<cudaReal> const & b,
                   const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                b.cArray()+beginIdB, n);
}

// Vector exponentiation, a[i] = exp(b[i]), kernel wrapper (cudaComplex).
__host__ void expV(DeviceArray<cudaComplex>& a, 
                   DeviceArray<cudaComplex> const & b,
                   const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                b.cArray()+beginIdB, n);
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

// Vector addition in-place, a[i] += b[i], GPU kernel (mixed).
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

// Vector addition in-place, a[i] += b, GPU kernel (mixed).
__global__ void _addEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x += b;
   }
}

// Vector addition in-place, a[i] += b[i], kernel wrapper (cudaReal).
__host__ void addEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector addition in-place, a[i] += b[i], kernel wrapper (cudaComplex).
__host__ void addEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector addition in-place, a[i] += b[i], kernel wrapper (mixed).
__host__ void addEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector addition in-place, a[i] += b, kernel wrapper (cudaReal).
__host__ void addEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector addition in-place, a[i] += b, kernel wrapper (cudaComplex).
__host__ void addEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector addition in-place, a[i] += b, kernel wrapper (mixed).
__host__ void addEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
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

// Vector subtraction in-place, a[i] -= b[i], GPU kernel (mixed).
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

// Vector subtraction in-place, a[i] -= b, GPU kernel (mixed).
__global__ void _subEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x -= b;
   }
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaReal).
__host__ void subEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (cudaComplex).
__host__ void subEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector subtraction in-place, a[i] -= b[i], kernel wrapper (mixed).
__host__ void subEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaReal).
__host__ void subEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector subtraction in-place, a[i] -= b, kernel wrapper (cudaComplex).
__host__ void subEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector subtraction in-place, a[i] -= b, kernel wrapper (mixed).
__host__ void subEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _subEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
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

// Vector multiplication in-place, a[i] *= b[i], GPU kernel (mixed).
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

// Vector multiplication in-place, a[i] *= b, GPU kernel (mixed).
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
__host__ void mulEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (cudaComplex).
__host__ void mulEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaComplex> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector multiplication in-place, a[i] *= b[i], kernel wrapper (mixed).
__host__ void
mulEqV(DeviceArray<cudaComplex>& a, DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaReal).
__host__ void mulEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (cudaComplex).
__host__ void mulEqS(DeviceArray<cudaComplex>& a, cudaComplex const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector multiplication in-place, a[i] *= b, kernel wrapper (mixed).
__host__ void mulEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _mulEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
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

// Vector division in-place, a[i] /= b[i], GPU kernel (mixed).
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

// Vector division in-place, a[i] /= b, GPU kernel (mixed).
__global__ void _divEqS(cudaComplex* a, cudaReal const b, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i].x /= b;
      a[i].y /= b;
   }
}

// Vector division in-place, a[i] /= b[i], kernel wrapper (cudaReal).
__host__ void divEqV(DeviceArray<cudaReal>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector division in-place, a[i] /= b[i], kernel wrapper (mixed).
__host__ void divEqV(DeviceArray<cudaComplex>& a, 
                     DeviceArray<cudaReal> const & b,
                     const int beginIdA, const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                  b.cArray()+beginIdB, n);
}

// Vector division in-place, a[i] /= b, kernel wrapper (cudaReal).
__host__ void divEqS(DeviceArray<cudaReal>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector division in-place, a[i] /= b, kernel wrapper (mixed).
__host__ void divEqS(DeviceArray<cudaComplex>& a, cudaReal const b,
                     const int beginIdA, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _divEqS<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, b, n);
}

// Vector addition in-place w/ coefficient, a[i] += b[i] * c, GPU kernel.
__global__ void _addEqVc(cudaReal* a, cudaReal const * b, 
                         cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] += b[i] * c;
   }
}

// Vector addition in-place w/ coefficient, a[i] += b[i] * c, kernel wrapper.
__host__ void addEqVc(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaReal> const & b, 
                      cudaReal const c, const int beginIdA, 
                      const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _addEqVc<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                   b.cArray()+beginIdB, c, n);
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
__host__ void sqNormV(DeviceArray<cudaReal>& a, 
                      DeviceArray<cudaComplex> const & b,
                      const int beginIdA, const int beginIdB, 
                      const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _sqNormV<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                   b.cArray()+beginIdB, n);
}

// Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), GPU kernel.
__global__ void _expVc(cudaReal* a, cudaReal const * b, 
                       cudaReal const c, const int n)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < n; i += nThreads) {
      a[i] = exp(b[i] * c);
   }
}

// Vector exponentiation w/ coefficient, a[i] = exp(b[i]*c), kernel wrapper.
__host__ void expVc(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b, 
                    cudaReal const c, const int beginIdA, 
                    const int beginIdB, const int n)
{
   UTIL_CHECK(a.capacity() >= n + beginIdA);
   UTIL_CHECK(b.capacity() >= n + beginIdB);
   
   // GPU resources
   int nBlocks, nThreads;
   ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

   // Launch kernel
   _expVc<<<nBlocks, nThreads>>>(a.cArray()+beginIdA, 
                                 b.cArray()+beginIdB, c, n);
}

// Add more than 2 vectors pointwise, GPU kernel.
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

// Add more than 2 vectors pointwise, kernel wrapper.
__host__ void addVMany(DeviceArray<cudaReal>& a, 
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
}

// Multiply more than 2 vectors pointwise, GPU kernel.
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

// Multiply more than 2 vectors pointwise, kernel wrapper.
__host__ void mulVMany(DeviceArray<cudaReal>& a, 
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
}

}
}
