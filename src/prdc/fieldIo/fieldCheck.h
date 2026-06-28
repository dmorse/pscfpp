#ifndef PRDC_FIELD_CHECK_H
#define PRDC_FIELD_CHECK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>   // Template with a default parameter

// Forward class declarations
namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf;

   /**
   * Allocate a DArray of fields.
   *
   * \param fields  DArray of fields of type FT
   * \param n  number of inner arrays
   * \param dimension  mesh dimensions
   */
   template <int D, class FT>
   void allocateFields(DArray<FT>& fields, int n, 
                       IntVec<D> const& dimension);

   /**
   * Check allocation of a single field, allocate if necessary.
   *
   * Template parameter FT is a field type, such as RField<D> or
   * RFieldDft<D>, that has an allocate function that takes an
   * IntVec<D> of mesh dimensions, and a dimensions function that
   * returns an IntVec<D> of mesh dimensions.
   *
   * On successful exit, the mesh dimensions for the field are
   * equal to those given in parameter meshDimensions.
   *
   * If the field is allocated on entry, the above condition is
   * checked, and an Exception is thrown if is violated.  If the
   * fields is not allocated on entry, it is allocated with the
   * correct dimensions.
   *
   * \ingroup Prdc_Field_Module
   *
   * \param field  field object of type FT
   * \param dimensions  required mesh dimensions
   */
   template <int D, class FT>
   void checkAllocateField(FT& field, IntVec<D> const& dimensions);

   /**
   * Check allocation of an array of fields, allocate if necessary.
   *
   * Template parameter FT is a field type, such as RField<D> or
   * RFieldDft<D>, that has an allocate function that takes an
   * IntVec<D> of mesh dimensions, and a dimensions function that
   * returns an IntVec<D> of mesh dimensions.
   *
   * On successful exit, the capacity of the DArray fields is equal
   * to nMonomer, and the mesh dimensions for each field are equal to
   * those given in parameter meshDimensions.
   *
   * If the fields array is allocated on entry, the above conditions
   * are checked, and an Exception is thrown if any are not met.
   * If the arrray is not allocated on entry, it is allocated with
   * the required dimensions.
   *
   * \ingroup Prdc_Field_Module
   *
   * \param fields  DArray of fields of type FT (in/out)
   * \param nMonomer  required number of monomer types (in)
   * \param dimensions  required mesh dimensions (in)
   */
   template <int D, class FT>
   void checkAllocateFields(DArray< FT >& fields,
                            int nMonomer,
                            IntVec<D> const& dimensions);

   /**
   * Inspect dimensions of a DArray of fields, each of type FT.
   *
   * Template parameter AT is an allocatable array type, such as
   * Util::DArray<double> or Pscf::HostDArray<double>, that has an allocate
   * function that takes an integer capacity as its only parameter.
   *
   * An Exception is thrown if fields is not allocated, or if the fields
   * do not all have the same positive list of dimensions.
   *
   * \ingroup Prdc_Field_Module
   *
   * \param fields  DArray of fields (in)
   * \param nMonomer  number of fields, or monomer types (out)
   * \param dimensions  dimensions mesh for each field (out)
   */
   template <int D, class FT>
   void inspectFields(DArray<FT> const & fields,
                      int & nMonomer,
                      IntVec<D> & dimensions);

   // Arrays of arrays 

   /**
   * Allocate an array of arrays.
   *
   * \param arrays  DArray of arrays
   * \param n  number of inner arrays
   * \param capacity  capacity of each inner array
   */
   template <class AT>
   void allocateArrays(DArray<AT>& arrays, int n, int capacity);

   /**
   * Check allocation of a DArray of 1D arrays, allocate if necessary.
   *
   * Template parameter AT is an allocatable array type, such as
   * Util::DArray<double> or Pscf::HostDArray<double>, that has an allocate
   * function that takes an integer capacity as its only parameter.
   *
   * On successful exit, the capacity of the arrays container is equal to
   * the parameter nMonomer, and the capacity for each element of type AT
   * is equal to the parameter capacity.
   *
   * If the container arrays is allocated on entry, the above conditions are
   * checked, and an Exception is thrown if any are violated.  If arrays is
   * not allocated on entry, it is allocated with the required dimensions.
   *
   * \ingroup Prdc_Field_Module
   *
   * \param arrays  DArray of arrays, each of type AT (in/out)
   * \param nMonomer  required number of arrays, or monomer types (in)
   * \param capacity  required capacity of each array (in)
   */
   template <int D, class AT>
   void checkAllocateArrays(DArray< AT >& arrays,
                            int nMonomer,
                            int capacity);

   /**
   * Inspect dimensions of a DArray of 1D arrays, each of type AT.
   *
   * An Exception is thrown if the arrays container is not allocated, or
   * if the 1D arrays do not all have the same positive capacity.
   *
   * \ingroup Prdc_Field_Module
   *
   * \param arrays  DArray of arrays (in)
   * \param nMonomer  number of arrays, or monomer types (out)
   * \param capacity  capacity of each array (out)
   */
   template <class AT>
   void inspectArrays(DArray<AT> const & arrays,
                      int & nMonomer,
                      int & capacity);

   /**
   * Copy a DArray of 1D arrays.
   *
   * The input and output containers may contain the same number of
   * inner arrays, each of the same positive capacity. 
   *
   * \ingroup Prdc_Field_Module
   *
   * \param out  DArray of arrays (lhs, out)
   * \param in  DArray of arrays (rhs, in)
   */
   template <class OAT, class IAT>
   void copyArrays(DArray<OAT>& out, DArray<IAT> const& in);

} // namespace Prdc
} // namespace Pscf
#include "fieldCheck.tpp"
#endif
