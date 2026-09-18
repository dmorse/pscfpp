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
   * Read-only psuedo-"host" array for use with C++ backend.
   *
   * This class template partial specialization may be used in template 
   * code in which data is copied from a psuedo-device array to a 
   * pseudo-host array when the device array is const and/or if the host 
   * array only requires read access. Specializations of this template
   * are directly derived from Util::ConstArrayView<Data>, and indirectly 
   * derived from Util::ConstArray<Data>.
   *
   * Construction or assignment (operator =)from a DeviceArray<Data,CPT> 
   * to a ConstHostArray<Data,CPT> creates a shallow read-only copy of a
   * C array that is owned by the device array. Assignment does nothing 
   * if such an association already exists.  Because the relevant 
   * constructor and assignment operator each take a const reference to
   * a DeviceArray<Data,CPT> as a parameter, they can be used in contexts 
   * in which the device array is declared const. This class template
   * hides the associate member functions of ConstArray<Data> by declaring 
   * these functions as private, so that construction and assignment are
   * the only allowed methods of creating a shallow copy.
   * 
   * The use of construction or an assignment (=) operator to create a 
   * shallow copy allows this class to be used in template code to 
   * imitate the syntax of an actual device-to-host data copy that 
   * would be performed by a ConstHostArray<Data, CUT> in a template
   * specialization that is designed to use a GPU.
   *
   * An association with a device array can be released by calling the
   * inherited dissociate function, or will be released upon destruction.
   * Such an association must be released by one of these two methods 
   * before the device array that owns the associated data is 
   * de-allocated or destroyed. A reference counting system detects and
   * reports de-allocation of a source array that is still referred to
   * by one or more other array view containers.
   * 
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class ConstHostArray<Data, CPT> : public ConstArrayView<Data>
   {

   public:

      using typename ConstArrayView<Data>::ValueType;

      // Default constructor
      ConstHostArray() = default;

      /**
      * Copy construction from a DeviceArray.
      *
      * Create an association (shallow copy) with memory owned by
      * the pre-existing device array.
      *
      * \param other  array container 
      */
      ConstHostArray(DeviceArray<Data,CPT> const & other);

      // Prohibit copy construction from another ConstHostArray.
      ConstHostArray(ConstHostArray<Data,CPT> const & other) = default;

      // Destructor
      ~ConstHostArray() = default;

      // Prohibit assignment from another ConstHostArray.
      ConstHostArray<Data,CPT>&
      operator = (ConstHostArray<Data,CPT> const&) = delete;

      /**
      * Create read-only association with a DeviceArray, if needed.
      *
      * If this is already associated with other (the device array), 
      * do nothing and return. Otherwise, create an association of this 
      * with the other device array (i.e., create a shallow copy).
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is associated with other with wrong size
      *
      * \param other  array container on RHS of assigment (input)
      */
      ConstHostArray<Data,CPT>& operator = (DeviceArray<Data,CPT> const & other);

   private:

      // Hide associate functions defined by base class.
      using ConstArrayView<Data>::associate;

   };

} // namespace Pscf

#include "DeviceArray.h"
namespace Pscf {

  /*
  * Copy construction from a DeviceArray.
  *
  * Creates an association with a DeviceArray, which must own the data.
  */
  template <typename Data>
  ConstHostArray<Data,CPT>::ConstHostArray(
                                 DeviceArray<Data,CPT> const & other)
   : ConstArrayView<Data>()
  {  ConstArrayView<Data>::associate(other); }

  /*
  * Assignment from a DeviceArray.
  *
  * Creates an association with a DeviceArray, which must own the data.
  */
  template <typename Data>
  ConstHostArray<Data,CPT>&
  ConstHostArray<Data,CPT>::operator = (DeviceArray<Data,CPT> const & other)
  {
     UTIL_CHECK(other.isAllocated());
     Data const * data = ConstArrayView<Data>::cArray();
     if ((bool)data && other.cArray() == data) {
        // If this is already associated with other, do nothing
        UTIL_CHECK(other.capacity() == ConstArrayView<Data>::size());
     } else {
        // Otherwise, attempt to create an association
        // Attempt fails if this is associated with a different array.
        ConstArrayView<Data>::associate(other);
     }
     return *this;
  }

}
#endif
