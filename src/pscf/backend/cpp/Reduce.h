#ifndef PSCF_CPU_REDUCE_H
#define PSCF_CPU_REDUCE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Functions that perform array reductions on a CPU.
   *
   * A reduction is any operation that involves reducing all of 
   * the elements of an array or set of arrays to a single scalar.  
   * Examples include taking the sum or finding the maximum of all 
   * array elements, or taking an inner product of two arrays.
   *
   * \defgroup Pscf_Cpu_Reduce_Module Reduce (CPU)
   * \ingroup Pscf_Cpu_Module
   */

   namespace Reduce {

      /*
      *  This file and Reduce.cpp contain declarations and definitions
      *  for reduction functions that operate real arrays.  Corresponding 
      *  declarations and definitions for reductions of complex-valued
      *  arrays are in files ReduceCx.h and ReduceCx.cpp, respectively. 
      */

      // Summation

      /**
      * Compute sum of array elements (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      * \return  sum of all array elements
      */
      double sum(Array<double> const & in);

      /**
      * Compute sum of elements of an array slice (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
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
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double sumSq(Array<double> const & in);

      // Summation of products

      /**
      * Compute Euclidean inner product of two real arrays .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param a  first input array
      * \param b  second input array
      */
      double innerProduct(Array<double> const & a,
                          Array<double> const & b);

      // Maxima 

      /**
      * Get maximum of array elements (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      */
      double max(Array<double> const & in);

      /**
      * Get value of maximum element in an array slice (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
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
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      * \return  maximum absolute value
      */
      double maxAbs(Array<double> const & in);

      // Minima

      /**
      * Get minimum of array elements .
      *
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      * \return  value of minimum element
      */
      double min(Array<double> const & in);

      /**
      * Get value of minimum element in an array slice (real).
      *
      * \ingroup Pscf_Cpu_Reduce_Module
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
      * \ingroup Pscf_Cpu_Reduce_Module
      *
      * \param in  input array
      * \return  value of minimum absolute value of all elements
      */
      double minAbs(Array<double> const & in);

   } // namespace Reduce
} // namespace Pscf
#endif
