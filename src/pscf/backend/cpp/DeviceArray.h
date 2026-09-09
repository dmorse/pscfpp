#ifndef PSCF_DEVICE_ARRAY_CP_H
#define PSCF_DEVICE_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FftwDRArray.h"                  // base class

// Forward declarations
namespace Pscf {
   template <typename Data, typename T> class DeviceArray;
   template <typename Data, typename T> class HostArray;
   class CPT;
}

namespace Pscf {

   using namespace Util;

   /**
   * Pseudo "device" array for use with C++ backend.
   *
   * Derived from FftwDRArray, and largely equivalent. 
   * 
   * The key difference from the base class is that assignment 
   * (operator = ) from a HostArray<Data,CTP> to a DeviceArray<Data,CTP> 
   * creates a shallow copy of the HostArray (a shared pointer) rather 
   * than a deep copy. This allows the creation of a shallow copy to be
   * used to imitate the syntax of an actual host-to-device data copy in 
   * template code that must work with either backend.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class DeviceArray<Data, CPT> : public FftwDRArray<Data>
   {

   public:

      using typename FftwDRArray<Data>::ValueType;

      // Default constructor
      DeviceArray() = default;

      // Copy constructor
      DeviceArray(DeviceArray<Data,CPT> const & other) = default;

      // Destructor
      ~DeviceArray() = default;

      // Assignment
      DeviceArray<Data,CPT>& 
      operator = (DeviceArray<Data,CPT> const&) = default;

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      *
      * \param capacity number of elements to allocate
      */
      DeviceArray(int capacity)
       : FftwDRArray<Data>(capacity)
      {}

      /**
      * Create association with a HostArray, or do nothing if associated.
      *
      * If both arrays are allocated and refer to the same memory block on 
      * entry, do nothing and return. Otherwise, if this is not allocated,
      * create an association of this with the HostArray (i.e., create
      * a shallow copy).
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is allocated and not associated with other
      *
      * \param other  array container on RHS of assigment (input)
      */
      DeviceArray<Data,CPT>& operator = (HostArray<Data,CPT> & other);

   };

} // namespace Pscf

#include "HostArray.h"
namespace Pscf {

  /*
  * Create an association with a HostArray, or do nothing if associated.
  */
  template <typename Data>
  DeviceArray<Data,CPT>& 
  DeviceArray<Data,CPT>::operator = (HostArray<Data,CPT> & other)
  {
     Data* data = Array<Data>::data_;
     bool same = (bool)data && data == other.cArray();
     if (same) {
        UTIL_CHECK(Array<Data>::capacity() == other.capacity());
     } else {
        FftwDRArray<Data>::associate(other); 
     }
     return *this;
  }

}
#endif
