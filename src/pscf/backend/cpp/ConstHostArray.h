#ifndef PSCF_CONST_HOST_ARRAY_CP_H
#define PSCF_CONST_HOST_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstArrayView.h>  // base class

// Forward declarations
namespace Pscf {
   template <typename Data, typename T> class ConstHostArray;
   template <typename Data, typename T> class DeviceArray;
   class CPT;
}

namespace Pscf {

   using namespace Util;

   /**
   * Read-only  psuedo-"host" array for use with C++ backend.
   *
   * This class is used for in template code in which data is copied from
   * a psuedo-device array to a pseudo-host array if the device array is
   * const and/or if the host array may be treated read only. It is
   * derived from Util::ArrayView, and is equivalent to this base class
   * except for the definition of a specialized assignment operator.
   *
   * Assignment from a DeviceArray<Data,CTP> const reference to a
   * ConstHostArray<Data,CTP> (operator =) creates a shallow read-only
   * copy of memory memory owned by the DeviceArray, as a Data const *
   * pointer owned by the ConstHostArray, or does nothing if such an
   * association already exists. This allows inexpensive creation of a
   * shallow read-only copy to be used to imitate the syntax of an
   * actual device-to-host data copy in template code that must work with
   * both CPU or GPU backends, without an unnecessary deep copy. Because
   * this assignment operator takes a const reference as its input, it
   * can can be used in contexts in which the device array is const.
   *
   * After an instance of this class creates an association with a
   * DeviceArray, this object must either be destroyed or explicitly
   * call ConstArrayView::dissociate before the associated DeviceArray
   * is de-allocated or destroyed. The destructor automatically releases
   * the association.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class ConstHostArray<Data, CPT> : public Util::ConstArrayView<Data>
   {

   public:

      using typename ConstArrayView<Data>::ValueType;

      // Default constructor
      ConstHostArray() = default;

      // Destructor
      ~ConstHostArray() = default;

      // Prohibit copy construction.
      ConstHostArray(ConstHostArray<Data,CPT> const & other) = default;

      // Prohibit assignment.
      ConstHostArray<Data,CPT>&
      operator = (ConstHostArray<Data,CPT> const&) = delete;

      /**
      * Create read-only association with a DeviceArray, if needed.
      *
      * If both arrays are allocated and point to the same memory block on
      * entry, do nothing and return. Otherwise, if this is not allocated,
      * create an association of this with the DeviceArray (i.e., create
      * a shallow copy).
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is allocated and not associated with other
      *
      * \param other  array container on RHS of assigment (input)
      */
      ConstHostArray<Data,CPT>& operator = (DeviceArray<Data,CPT> & other);

   };

} // namespace Pscf

#include "DeviceArray.h"
namespace Pscf {

  /*
  * Create an association with a DeviceArray, or do nothing if associated.
  */
  template <typename Data>
  ConstHostArray<Data,CPT>&
  ConstHostArray<Data,CPT>::operator = (DeviceArray<Data,CPT> & other)
  {
     const Data* data = ConstArrayView<Data>::cArray();
     if ((bool)data && other.cArray() == data) {
        UTIL_CHECK(other.capacity() == ConstArrayView<Data>::size());
        // If this is already associated with the other array, do nothing
     } else {
        // Otherwise, create an association
        ConstArrayView<Data>::associate(other);
     }
     return *this;
  }

}
#endif
