/*
* PSCF Package
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Reduce.h"
#include <util/containers/Array.h>
#include <cmath>

namespace Pscf {
namespace Prdc {
namespace Cuda {
namespace Reduce {

   /*
   * Get maximum of array elements.
   */
   double max(Array<double> const & in)
   {
      int n = in.capacity();
      double max = in[0];
      for (int i = 1; i < n; i++) {
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
      double min = in[0];
      for (int i = 1; i < n; i++) {
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
      double val;
      double min = std::abs(in[0]);
      for (int i = 1; i < n; i++) {
         val = std::abs(in[i]);
         if (val < min) min = val;
      }
      return min;
   }

   /*
   * Compute sum of array elements (GPU kernel wrapper).
   */
   double sum(Array<double> const & in)
   {
      int n = in.capacity();
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += in[i];
      }
      return sum;
   }

   /*
   * Compute inner product of two real arrays.
   */
   double innerProduct(Array<double> const & a,
                       Array<double> const & b)
   {
      int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
         sum += a[i]*b[i];
      }
      return sum;
   }

} // Reduce
} // Cpu
} // Prdc
} // Pscf
