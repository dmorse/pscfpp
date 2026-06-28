#ifndef PRDC_FIELD_CHECK_TPP
#define PRDC_FIELD_CHECK_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {

   // Arrays of fields 

   /*
   * Allocate an array of n fields of type FT.
   */
   template <int D, class FT>
   void allocateFields(DArray<FT>& fields, 
                       int n, 
                       IntVec<D> const& dimension)
   {
      UTIL_CHECK(!fields.isAllocated());
      fields.allocate(n);
      for (int i = 0; i < n; ++i) {
         fields[i].allocate(dimension);	 
      }
   }

   /*
   * Check allocation of a single field of type FT, allocate if needed.
   */
   template <int D, class FT>
   void checkAllocateField(FT& field, 
                           IntVec<D> const& dimensions)
   {
      if (field.isAllocated()) {
         UTIL_CHECK(field.meshDimensions() == dimensions);
      } else {
         field.allocate(dimensions);
      }
   }

   /*
   * Check an array of n fields of type FT, allocate if needed.
   */
   template <int D, class FT>
   void checkAllocateFields(DArray<FT>& fields,
                            int nMonomer, 
                            IntVec<D> const& dimensions)
   {
      if (fields.isAllocated()) {
         int nMonomerFields = fields.capacity();
         UTIL_CHECK(nMonomerFields > 0)
         UTIL_CHECK(nMonomerFields == nMonomer)
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(fields[i].isAllocated());
            UTIL_CHECK(fields[i].meshDimensions() == dimensions);
         }
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(dimensions);
         }
      }
   }

   /*
   * Inspect dimensions of a DArray of n fields of type FT.
   */
   template <int D, class FT>
   void inspectFields(DArray<FT> const& fields,
                      int & nMonomer,
                      IntVec<D> & dimensions)
   {
      UTIL_CHECK(fields.isAllocated());
      nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      dimensions = fields[0].meshDimensions();
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].isAllocated());
         UTIL_CHECK(fields[i].meshDimensions() == dimensions);
      }
   }

   // Array of arrays 

   template <class AT>
   void allocateArrays(DArray<AT>& arrays, int n, int capacity)
   {
      UTIL_CHECK(!arrays.isAllocated());
      arrays.allocate(n);
      for (int i = 0; i < n; ++i) {
         arrays[i].allocate(capacity);	 
      }
   }

   /*
   * Check allocation of an array of arrays of type AT, allocate if needed.
   */
   template <class AT>
   void checkAllocateArrays(DArray<AT>& arrays,
                            int nMonomer, 
                            int capacity)
   {
      if (arrays.isAllocated()) {
         int nMonomerArrays = arrays.capacity();
         UTIL_CHECK(nMonomerArrays > 0)
         UTIL_CHECK(nMonomerArrays == nMonomer)
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(arrays[i].isAllocated());
            UTIL_CHECK(arrays[i].capacity() == capacity);
         }
      } else {
         arrays.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            arrays[i].allocate(capacity);
         }
      }
   }

   /*
   * Inspect dimensions of an array of arrays of type AT.
   */
   template <class AT>
   void inspectArrays(DArray<AT> const& arrays,
                    int & nMonomer,
                    int & capacity)
   {
      UTIL_CHECK(arrays.isAllocated());
      nMonomer = arrays.capacity();
      UTIL_CHECK(nMonomer > 0);

      capacity = arrays[0].capacity();
      UTIL_CHECK(capacity > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(arrays[i].isAllocated());
         UTIL_CHECK(arrays[i].capacity() == capacity);
      }
   }

   template <class OAT, class IAT>
   void copyArrays(DArray<OAT>& out, DArray<IAT> const& in) 
   {
      int n = in.capacity();
      UTIL_CHECK(out.capacity() == n);
      for (int i = 0; i < n; ++i) {
	 UTIL_CHECK(in[i].isAllocated());
	 UTIL_CHECK(out[i].isAllocated());
	 UTIL_CHECK(in[i].capacity() == out[i].capacity());
         out[i] = in[i];
      }
   }

} // namespace Prdc
} // namespace Pscf
#endif
