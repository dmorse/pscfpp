#ifndef PSCF_CUDA_REDUCE_H
#define PSCF_CUDA_REDUCE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>
#include <complex>

namespace Pscf {

   /**
   * Functions that perform parallel reductions on a GPU.
   *
   * A reduction is any operation that involves reducing all of the
   * elements of one or more arrays to a single scalar result. Examples
   * include taking the computation of the sum or the maximum of all
   * elements of a single array, or of a Euclidean inner product of
   * two real arrays.
   *
   * A public function is provided for each reduction operation, which
   * takes one or more DeviceArray containers as inputs and returns a
   * single value. Each wrapper function is a standard C++ function that
   * may be called on the host CPU. The implementation of each such
   * wrapper function calls one of the device-wide algorithms defined
   * in the CUB library that is distributed as part of the CUDA Core
   * Compute library collection provided with the CUDA development kit.
   *
   * Reduction functions that return a real number return a value of
   * type cudaReal. Functions that return a complex number instead return
   * a value of type std::complex<cudaReal>. Complex values are returned
   * as std::complex<cudaReal> rather than cudaComplex to allow similar
   * interfaces to be used by reduction functions declared here and by
   * analogous functions defined in src/pscf/cpu/Reduced.h that take
   * take standard Array containers as inputs and perform calculations
   * on the CPU.
   *
   * \defgroup Pscf_Cuda_Reduce_Module Reduce (GPU)
   * \ingroup Pscf_Cuda_Module
   */

   /**
   * Reduction operations performed on a CPU or GPU.
   *
   * \ingroup Pscf_Math_Module
   */
   namespace Reduce {

      // Summation

      /**
      * Return sum of elements of a real array.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of elements
      */
      cudaReal sum(DeviceArray<cudaReal> const & in);

      /**
      * Return sum of elements of a real array slice.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \param begin  index of first element of slice
      * \param end  index one past last element
      * \return  sum of elements
      */
      cudaReal sum(DeviceArray<cudaReal> const & in, int begin, int end);

      /**
      * Return sum of all elements of a complex array.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  complex input array
      * \return  complex sum of elements
      */
      std::complex<cudaReal> sum(DeviceArray<cudaComplex> const & in);

      /**
      * Return sum of elements of a complex array slice.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \param begin  index of first element of slice
      * \param end  index one past last element of slice
      * \return  sum of elements
      */
      std::complex<cudaReal> sum(DeviceArray<cudaComplex> const & in,
                                 int begin, int end);

      // Sum of squares and array products

      /**
      * Return sum of squares of elements of a real array.
      *
      * This function returns the square of the Euclidean norm for a
      * real array.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of squares of elements
      */
      cudaReal sumSq(DeviceArray<cudaReal> const & in);

      /**
      * Return sum of squares of elements of a complex array.
      *
      * This function returns the complex sum of complex squares of
      * elements. This is not the square of the Hilbert space norm.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of complex squares of elements
      */
      std::complex<cudaReal> sumSq(DeviceArray<cudaComplex> const & in);

      /**
      * Return sum of squared magnitudes of elements of a complex array.
      *
      * This function returns the real sum of the squared absolute
      * magnitude of all elements of a complex array. This is the square
      * of the conventional Hilbert space norm.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  sum of squared absolute magnitude of elements
      */
      cudaReal sumSqAbs(DeviceArray<cudaComplex> const & in);

      /**
      * Return the inner product of two real arrays.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param a  first real input array
      * \param b  second real input array
      * \return  Euclidean inner product
      */
      cudaReal innerProduct(DeviceArray<cudaReal> const & a,
                            DeviceArray<cudaReal> const & b);

      // Maxima - real array inputs

      /**
      * Return maximum of all real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  value of maximum element
      */
      cudaReal max(DeviceArray<cudaReal> const & in);

      /**
      * Return maximum of elements of a real array slice.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \param begin  index of first element of slice
      * \param end  index one past last element
      * \return  value of maximum element
      */
      cudaReal max(DeviceArray<cudaReal> const & in, int begin, int end);

      /**
      * Get maximum absolute magnitude of real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  magnitude of element of maximum magnitude
      */
      cudaReal maxAbs(DeviceArray<cudaReal> const & in);

      // Minima - real array inputs

      /**
      * Return minimum of all real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  value of minimum element
      */
      cudaReal min(DeviceArray<cudaReal> const & in);

      /**
      * Return minimum of elements of a real array slice.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \param begin  index of first element of slice
      * \param end  index one past last element
      * \return  value of minimum element
      */
      cudaReal min(DeviceArray<cudaReal> const & in, int begin, int end);

      /**
      * Return minimum absolute magnitude of real array elements.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      *
      * \param in  real input array
      * \return  magnitude of element of minimum magnitude
      */
      cudaReal minAbs(DeviceArray<cudaReal> const & in);

      // Memory management

      /**
      * Free any private work space currently allocated for reductions.
      *
      * \ingroup Pscf_Cuda_Reduce_Module
      */
      void freeWorkSpace();

   } // namespace Reduce
} // namespace Pscf
#endif
