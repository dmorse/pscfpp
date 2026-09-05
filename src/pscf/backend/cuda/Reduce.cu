/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Reduce.h"
#include "VecOp.h"
#include <pscf/backend/cuda/ThreadArray.h>
#include <pscf/backend/cuda/HostDArray.h>
#include <pscf/backend/cuda/DeviceMemory.h>
#include <pscf/backend/cuda/cudaErrorCheck.h>
#include <pscf/backend/cuda/complex.h>
#include <util/misc/Log.h>

//#include <thrust/reduce.h>
//#include <thrust/device_vector.h>
//#include <cub/device/device_reduce.cuh>
#include <cub/cub.cuh>

#include <complex>
#include <cmath>

namespace Pscf {
namespace Reduce {

   /*
   * Reduction functions defined here use device-wide algorithms defined
   * in the CUB library. The CUB library is part of the CUDA Core Compute
   * Library collection that is distributed with the CUDA development kit.
   *
   * The CUB reduction functions require the user to provide workspace.
   * To minimize time spent repeatedly allocating this workspace, the
   * reduceSpace static variable is used as a memory block that grows as
   * needed, but that can be re-used for operations that require a work
   * space of size less than or equal to that of the largest previous
   * required work space.
   *
   * The static variable transformSpace is used as a memory block to
   * store the results of vector operations (such as taking the absolute
   * magnitude or square of all elements) that are sometimes applied
   * before applying a reduction operation.
   *
   * Both memory pools may be de-allocated by Reduce::freeWorkSpace(),
   * but otherwise persist in memory for use as needed.
   */

   // Memory used as workspace by CUB reduction functions
   static DeviceMemory reduceSpace_{};

   // Memory used as workspace for vector operations (e.g, max or square)
   static DeviceMemory transformSpace_{};

   // Anonymous namespace, for entities that are only used in this file
   namespace {

      // Complex addition functor.
      struct addComplexFunctor {

         __host__ __device__ inline
         cudaComplex operator() (cudaComplex const & a,
                                 cudaComplex const & b) const
         {
            cudaComplex result;
            result.x = a.x + b.x;
            result.y = a.y + b.y;
            return result;
         }

      };

      /*
      * Compute sum of elements of a real array.
      */
      cudaReal sum(cudaReal* inPtr, const int n)
      {
         UTIL_CHECK(n > 0);

         // Create output pointer
         DeviceArray<cudaReal> out(1);
         cudaReal* outPtr = out.cArray();

         // Determine size of required workspace, allocate if needed
         size_t workSize = 0;
         cudaError_t error;
         error = cub::DeviceReduce::Sum(nullptr, workSize,
                                        inPtr, outPtr, n);
         UTIL_CHECK(error == cudaSuccess);
         reduceSpace_.resize(workSize);
         UTIL_CHECK(reduceSpace_.capacity() >= workSize);

         // Perform reduction
         error = cub::DeviceReduce::Sum(reduceSpace_.cArray(), workSize,
                                        inPtr, outPtr, n);
         UTIL_CHECK(error == cudaSuccess);

         // Copy to host and return value
         HostDArray<cudaReal> out_h;
         out_h.allocate(1);
         out_h = out;
         return out_h[0];
      }

      /*
      * Compute sum of elements of a complex array.
      *
      * Implementation uses Nvidia CUB Library functions.
      */
      std::complex<cudaReal> sum(cudaComplex* inPtr, const int n)
      {
         UTIL_CHECK(n > 0);

         // Create output pointer
         DeviceArray<cudaComplex> out(1);
         cudaComplex* outPtr = out.cArray();

         // Determine size of required workspace and allocate if necessary
         size_t workSize = 0;
         cudaError_t error;
         auto op = addComplexFunctor{};
         cudaComplex init = makeComplex(0.0, 0.0);
         error = cub::DeviceReduce::Reduce(nullptr, workSize,
                                           inPtr, outPtr, n, op, init);
         UTIL_CHECK(error == cudaSuccess);
         reduceSpace_.resize(workSize);

         // Perform reduction
         error = cub::DeviceReduce::Reduce(reduceSpace_.cArray(), workSize,
                                           inPtr, outPtr, n, op, init);
         UTIL_CHECK(error == cudaSuccess);

         // Copy to host and return value
         HostDArray<cudaComplex> out_h(1);
         out_h = out;
         return std::complex<cudaReal>(out_h[0].x, out_h[0].y);
      }

      /*
      * Compute max of elements of a real array.
      */
      cudaReal max(cudaReal* inPtr, const int n)
      {
         UTIL_CHECK(n > 0);

         // Create input and output pointers
         DeviceArray<cudaReal> out(1);
         cudaReal* outPtr = out.cArray();

         // Determine size of required workspace, allocated if needed
         size_t workSize = 0;
         cudaError_t error;
         error = cub::DeviceReduce::Max(nullptr, workSize,
                                        inPtr, outPtr, n);
         UTIL_CHECK(error == cudaSuccess);
         reduceSpace_.resize(workSize);

         // Perform reduction
         error = cub::DeviceReduce::Max(reduceSpace_.cArray(), workSize,
                                        inPtr, outPtr, n);
         UTIL_CHECK(error == cudaSuccess);

         // Copy to host and return value
         HostDArray<cudaReal> out_h(1);
         out_h = out;
         return out_h[0];
      }

      /*
      * Compute min of elements of a real array.
      */
      cudaReal min(cudaReal* inPtr, const int n)
      {
         UTIL_CHECK(n > 0);

         // Create output pointer
         DeviceArray<cudaReal> out(1);
         cudaReal* outPtr = out.cArray();

         // Determine size of required workspace, allocated if needed
         size_t workSize = 0;
         cudaError_t error;
         error = cub::DeviceReduce::Min(nullptr, workSize,
                                        inPtr, outPtr, n);
         UTIL_CHECK(error == cudaSuccess);
         reduceSpace_.resize(workSize);

         // Perform reduction
         error = cub::DeviceReduce::Min(reduceSpace_.cArray(), workSize,
                                        inPtr, outPtr, n);
         UTIL_CHECK(error == cudaSuccess);

         // Copy to host and return value
         HostDArray<cudaReal> out_h(1);
         out_h = out;
         return out_h[0];
      }

   } // anonmyous namespace Pscf::Reduce::(unnamed)

   // Memory management

   void freeWorkSpace()
   {
      reduceSpace_.deallocate();
      transformSpace_.deallocate();
   }

   // Reduction function definitions

   /*
   * Compute sum of elements of a real array.
   *
   * This implementation uses Nvidia CUB Library functions.
   */
   cudaReal sum(DeviceArray<cudaReal> const & a)
   {
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      return sum(inPtr, n);
   }

   /*
   * Compute sum of elements of a real array slice.
   *
   * This implementation uses Nvidia CUB Library functions.
   */
   cudaReal sum(DeviceArray<cudaReal> const & a, int begin, int end)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= a.capacity());
      const int n = end - begin;
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      inPtr += begin;
      return sum(inPtr, n);
   }

   /*
   * Compute sum of elements of a complex array.
   *
   * Implementation uses Nvidia CUB Library functions.
   */
   std::complex<cudaReal> sum(DeviceArray<cudaComplex> const & a)
   {
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      cudaComplex* inPtr = const_cast<cudaComplex*>( a.cArray() );
      return sum(inPtr, n);
   }

   /*
   * Compute sum of elements of a complex array slice.
   *
   * This implementation uses Nvidia CUB Library functions.
   */
   std::complex<cudaReal> sum(DeviceArray<cudaComplex> const & a,
                              int begin, int end)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= a.capacity());
      const int n = end - begin;
      UTIL_CHECK(n > 0);

      cudaComplex* inPtr = const_cast<cudaComplex*>( a.cArray() );
      inPtr += begin;
      return sum(inPtr, n);
   }

   /*
   * Return the sum of squares of elements of a real array.
   */
   cudaReal sumSq(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of element-wise squares
      VecOp::sqV(temp, in);

      cudaReal result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }

   /*
   * Return the sum of squares of elements of a complex array.
   */
   std::complex<cudaReal> sumSq(DeviceArray<cudaComplex> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaComplex);
      transformSpace_.resize(workSize);
      DeviceArray<cudaComplex> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of element-wise squares
      VecOp::sqV(temp, in);

      // Compute sum of array temp of element-wise squares
      std::complex<cudaReal> result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }

   /*
   * Return the sum of square absolute magnitudes for a complex array.
   */
   cudaReal sumSqAbs(DeviceArray<cudaComplex> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of element-wise square absolute magnitudes
      VecOp::sqAbsV(temp, in);

      // Compute sum of array temp of element-wise square magnitudes
      cudaReal result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }

   /*
   * Compute inner product of two real arrays.
   */
   cudaReal innerProduct(DeviceArray<cudaReal> const & a,
                         DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(b.isAllocated());
      UTIL_CHECK(a.capacity() == b.capacity());
      int n = a.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Perform element-wise multiplication v[i] = a[i] * b[i]
      VecOp::mulVV(temp, a, b);

      // Compute sum of element-wise products
      cudaReal result = Reduce::sum(temp);
      temp.dissociate();

      return result;
   }

   // Maxima

   /*
   * Compute max of elements of a real array.
   */
   cudaReal max(DeviceArray<cudaReal> const & a)
   {
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      return max(inPtr, n);
   }

   /*
   * Compute max of elements of a real array slice.
   */
   cudaReal max(DeviceArray<cudaReal> const & a, int begin, int end)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= a.capacity());
      const int n = end - begin;
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      inPtr += begin;
      return max(inPtr, n);
   }

   /*
   * Return the maximum absolute magnitude of array elements.
   */
   cudaReal maxAbs(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of absolute magnitudes
      VecOp::absV(temp, in);

      // Compute maximum of array temp of absolute magnitudes
      cudaReal result = Reduce::max(temp);
      temp.dissociate();

      return result;
   }

   // Minima

   /*
   * Compute minimum of all elements of a real array.
   */
   cudaReal min(DeviceArray<cudaReal> const & a)
   {
      UTIL_CHECK(a.isAllocated());
      const int n = a.capacity();
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      return min(inPtr, n);
   }

   /*
   * Compute minimum of elements of a real array slice.
   */
   cudaReal min(DeviceArray<cudaReal> const & a, int begin, int end)
   {
      UTIL_CHECK(a.isAllocated());
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= a.capacity());
      const int n = end - begin;
      UTIL_CHECK(n > 0);

      cudaReal* inPtr = const_cast<cudaReal*>( a.cArray() );
      inPtr += begin;
      return min(inPtr, n);
   }

   /*
   * Return the minimum absolute magnitude of array elements.
   */
   cudaReal minAbs(DeviceArray<cudaReal> const & in)
   {
      UTIL_CHECK(in.isAllocated());
      int n = in.capacity();

      // Set up temporary array for result of vector operation
      int workSize = n * sizeof(cudaReal);
      transformSpace_.resize(workSize);
      DeviceArray<cudaReal> temp;
      temp.associate(transformSpace_, n);

      // Compute an array of absolute magnitudes
      VecOp::absV(temp, in);

      // Compute minimum of array temp of absolute magnitudes
      cudaReal result = Reduce::min(temp);
      temp.dissociate();

      return result;
   }

} // namespace Reduce
} // namespace Pscf
