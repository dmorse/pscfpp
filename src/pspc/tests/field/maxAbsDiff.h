#ifndef PSPC_MAX_ABS_DIFF_H
#define PSPC_MAX_ABS_DIFF_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/container/DArray.h>


template <class ArrayT, typename absT = double>
absT maxAbsDiff(ArrayT const& a, ArrayT const& b)
{
   UTIL_CHECK(a.capacity() > 0);
   UTIL_CHECK(a.capacity() == b.capacity());
   int n = a.capacity();
   absT x;
   absT max = abs(a[0]);
   for (int i = 1; i < n; ++i) {
      x = abs(a[i] - b[i]); 
      if (x > max) max = x;
   }
   return max;
}

template <class ArrayT, typename absT = double>
absT maxAbsDiff(DArray<ArrayT> const& a, DArray<ArrayT> const& b)
{
   UTIL_CHECK(a.capacity() > 0);
   UTIL_CHECK(a.capacity() == b.capacity());
   UTIL_CHECK(a[0].capacity() > 0);
   int m = a.capacity();
   absT x;
   absT max = abs(a[0][0]);
   int i, j;
   for (i = 1; i < m; ++i) {
      n = a[i].capacity()
      UTIL_CHECK(n > 0);
      UTIL_CHECK(n == b[i].capacity());
      for (j = 0; j < n; ++j) {
         x = abs(a[i] - b[i]); 
         if (x > max) max = x;
      }
   }
   return max;
}
#endif
