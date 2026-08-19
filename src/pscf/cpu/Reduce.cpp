/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Reduce.h"
#include <util/containers/Array.h>
#include <cmath>

namespace Pscf {
namespace Reduce {

   // Summation

   /*
   * Compute the sum of array elements (real).
   */
   double sum(Array<double> const & in)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);

      // Kahan summation - to reduce accumulation of error
      double sum = 0.0;
      double err = 0.0;
      double tempVal, tempSum;
      for (int i = 0; i < n; ++i) {
         tempVal = in[i] - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;
      }

      #if 0
      // Simple summation
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += in[i];
      }
      #endif

      return sum;
   }

   /*
   * Compute the sum of array elements (real).
   */
   double sum(Array<double> const & in, int begin, int end)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= n);

      // Kahan summation - to reduce accumulation of error
      double sum = 0.0;
      double err = 0.0;
      double tempVal, tempSum;
      for (int i = begin; i < end; ++i) {
         tempVal = in[i] - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;
      }

      #if 0
      // Simple summation
      double sum = 0.0;
      for (int i = begin; i < end; i++) {
         sum += in[i];
      }
      #endif

      return sum;
   }

   /*
   * Compute the sum of squares of all array elements (real).
   */
   double sumSq(Array<double> const & in)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);

      // Kahan summation - to reduce accumulation of error
      double sum = 0.0;
      double err = 0.0;
      double x, tempVal, tempSum;
      for (int i = 0; i < n; ++i) {
         x = in[i];
         tempVal = x*x - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;
      }

      #if 0
      // Simple summation
      double val;
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         val = in[i];
         sum += val*val;
      }
      #endif

      return sum;
   }

   // Inner product

   /*
   * Compute Euclidean inner product of two arrays (real).
   */
   double innerProduct(Array<double> const & a,
                       Array<double> const & b)
   {
      int n = a.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(b.capacity() == n);

      // Kahan summation - to reduce accumulation of error
      double sum = 0.0;
      double err = 0.0;
      double tempVal, tempSum;
      for (int i = 0; i < n; ++i) {
         tempVal = a[i]*b[i] - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;
      }

      #if 0
      // Simple summation
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += a[i]*b[i];
      }
      #endif

      return sum;
   }

   // Maxima and minima

   /*
   * Get maximum of array elements.
   */
   double max(Array<double> const & in)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      double max = in[0];
      for (int i = 1; i < n; i++) {
         if (in[i] > max) max = in[i];
      }
      return max;
   }

   /*
   * Get maximum element of array slice.
   */
   double max(Array<double> const & in, int begin, int end)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= n);
      UTIL_CHECK(end > begin);
      double max = in[begin];
      for (int i = begin + 1; i < end; i++) {
         if (in[i] > max) max = in[i];
      }
      return max;
   }

   /*
   * Get maximum absolute magnitude of array elements.
   */
   double maxAbs(Array<double> const & in)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      double val;
      double max = std::abs(in[0]);
      for (int i = 1; i < n; i++) {
         val = std::abs(in[i]);
         if (val > max) max = val;
      }
      return max;
   }

   /*
   * Get minimum of array elements.
   */
   double min(Array<double> const & in)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      double min = in[0];
      for (int i = 1; i < n; i++) {
         if (in[i] < min) min = in[i];
      }
      return min;
   }

   /*
   * Get value of minimum element in an array slice.
   */
   double min(Array<double> const & in, int begin, int end)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      UTIL_CHECK(begin >= 0);
      UTIL_CHECK(end <= n);
      UTIL_CHECK(end > begin);
      double min = in[begin];
      for (int i = begin + 1; i < end; i++) {
         if (in[i] < min) min = in[i];
      }
      return min;
   }

   /*
   * Get minimum absolute magnitude of array elements.
   */
   double minAbs(Array<double> const & in)
   {
      int n = in.capacity();
      UTIL_CHECK(n > 0);
      double val;
      double min = std::abs(in[0]);
      for (int i = 1; i < n; i++) {
         val = std::abs(in[i]);
         if (val < min) min = val;
      }
      return min;
   }

} // Reduce
} // Pscf
