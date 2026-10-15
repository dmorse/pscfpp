#ifndef PSCF_CPU_REDUCE_H
#define PSCF_CPU_REDUCE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/complex.h>
#include <fftw3.h>

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Functions that perform reductions (C++ CPU backend).
   *
   * A reduction is a operation that takes one or more arrays as inputs 
   * and returns a scalar result, such as a sum of elements, a maximum 
   * or minimum element, or an inner product. 
   *
   * Complex reduction results are returned as std::complex<double>
   * values because the fftw_complex type is an alias for an array types,
   * double[2], that cannot be returned by value.
   *
   * \defgroup Pscf_Cpp_Reduce_Module Reduce (C++ backend)
   * \ingroup Pscf_Backend_Cpp_Module
   */

   /**
   * Functions that perform array reductions (C++ backend).
   */
   namespace Reduce {

      // Summation

      /**
      * Compute sum of array elements (real).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \return  sum of all array elements
      */
      double sum(Array<double> const & in);

      /**
      * Compute sum of elements of an array slice (real).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \param begin  index of first element in slice
      * \param end  index one past the last in slice
      * \return  sum of elements with indices [begin, end-1]
      */
      double sum(Array<double> const & in, int begin, int end);

      /**
      * Compute sum of of squares of array elements (real).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      */
      double sumSq(Array<double> const & in);

      /**
      * Compute sum of all elements of an array (complex).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param a  input array
      * \return sum of all elements of complex array a
      */
      std::complex<double> sum(Array<fftw_complex> const & a);

      /**
      * Compute sum of elements of an array slice (complex).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param a  input array
      * \param begin  index of first element of slice
      * \param end  index one past the last element
      * \return sum of elements with indices [begin, end-1]
      */
      std::complex<double> sum(Array<fftw_complex> const & a,
                               int begin, int end);

      /**
      * Compute sum of squares of elements of a complex array.
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param a  input array
      * \return complex sum of squares of elements of array a
      */
      std::complex<double> sumSq(Array<fftw_complex> const & a);

      // Summation of products

      /**
      * Compute Euclidean inner product of two real arrays .
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param a  first input array
      * \param b  second input array
      */
      double innerProduct(Array<double> const & a,
                          Array<double> const & b);

      /**
      * Compute sum of complex products of elements of two arrays (complex).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param a  first input array
      * \param b  second input array
      * \return sum of elementwise complex products a[i]*b[i]
      */
      std::complex<double> sumProduct(Array<fftw_complex> const & a,
                                      Array<fftw_complex> const & b);

      // Maxima 

      /**
      * Get maximum of array elements (real).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      */
      double max(Array<double> const & in);

      /**
      * Get value of maximum element in an array slice (real).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \param begin  index of first element in slice
      * \param end  index one past the last in slice
      * \return  maximum of elements with indices [begin, end-1]
      */
      double max(Array<double> const & in, int begin, int end);

      /**
      * Get maximum absolute magnitude of array elements .
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \return  maximum absolute value
      */
      double maxAbs(Array<double> const & in);

      // Minima

      /**
      * Get minimum of array elements .
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \return  value of minimum element
      */
      double min(Array<double> const & in);

      /**
      * Get value of minimum element in an array slice (real).
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \param begin  index of first element in slice
      * \param end  index one past the last in slice
      * \return  minimum of elements with indices [begin, end-1]
      */
      double min(Array<double> const & in, int begin, int end);

      /**
      * Get minimum absolute magnitude of array elements .
      *
      * \ingroup Pscf_Cpp_Reduce_Module
      *
      * \param in  input array
      * \return  value of minimum absolute value of all elements
      */
      double minAbs(Array<double> const & in);

   } // namespace Reduce
} // namespace Pscf
#endif
