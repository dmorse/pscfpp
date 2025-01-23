#ifndef PRDC_FIELD_ARRAY_UTIL_H
#define PRDC_FIELD_ARRAY_UTIL_H

/*
* PSCF Package 
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/types.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   template <class AT>
   void allocateArrays(DArray<AT>& arrays, int n, int capacity)
   {
      UTIL_CHECK(!arrays.isAllocated());
      arrays.allocate(n);
      for (int i = 0; i < n; ++i) {
         arrays[i].allocate(capacity);	 
      }
   }

   template <int D, class FT>
   void allocateFields(DArray<FT>& fields, int n, IntVec<D> const& dimension)
   {
      UTIL_CHECK(!fields.isAllocated());
      fields.allocate(n);
      for (int i = 0; i < n; ++i) {
         fields[i].allocate(dimension);	 
      }
   }

   template <class OAT, class IAT>
   void copyArrays(DArray<OAT>& out, DArray<IAT> const& in) 
   {
      int n = in.capacity();
      UTIL_CHECK(out.capacity() == n);
      for (int i = 0; i < n; ++i) {
	 UTIL_CHECK(in[i].capacity() == out[i].capacity());
         out[i] = in[i];
      }
   }

}
#endif
