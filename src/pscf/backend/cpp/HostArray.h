#ifndef PSCF_HOST_ARRAY_CP_H
#define PSCF_HOST_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FftwDRArray.h"                  // base class

// Forward declarations
namespace Pscf {
   template <typename Data, typename T> class HostArray;
   template <typename Data, typename T> class DeviceArray;
   class CPT;
}

namespace Pscf {

   using namespace Util;

   /**
   * Pseudo "host" array for use with C++ backend.
   *
   * Derived from FftwDRArray, and largely equivalent.
   *
   * The key difference from the FftwDRArray base class is that assignment
   * (operator = ) from a DeviceArray<Data,CTP> to a HostArray<Data,CTP>
   * creates a shallow copy of the DeviceArray (a shared pointer) rather
   * than a deep copy. This allows the inexpensive creation of a shallow 
   * copy to be used to imitate the syntax of an actual device-to-host data 
   * copy in template code that must work with CPU or GPU backend, without
   * the cost of an unnecessary copy.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class HostArray<Data, CPT> : public FftwDRArray<Data>
   {

   public:

      using typename FftwDRArray<Data>::ValueType;

      // Default constructor
      HostArray() = default;

      // Copy constructor
      HostArray(HostArray<Data,CPT> const & other) = default;

      // Destructor
      ~HostArray() = default;

      // Assignment
      HostArray<Data,CPT>& 
      operator = (HostArray<Data,CPT> const&) = default;

      /**
      * Allocating constructor.
      *
      * \param capacity number of elements to allocate
      */
      HostArray(int capacity)
       : FftwDRArray<Data>(capacity)
      {}

      /**
      * Create association with a DeviceArray, or do nothing if associated.
      *
      * If both arrays are allocated and point to the same address on 
      * entry, do nothing and return. Otherwise, if this is not allocated,
      * create an association of this with the DeviceArray (i.e., create
      * a shallow copy).
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is allocated and not associated with other
      *
      * \param other  array container on RHS of assigment (input)
      */
      HostArray<Data,CPT>& operator = (DeviceArray<Data,CPT> & other);

   };

} // namespace Pscf

#include "DeviceArray.h"
namespace Pscf {

  /*
  * Create an association with a DeviceArray, or do nothing if associated.
  */
  template <typename Data>
  HostArray<Data,CPT>&
  HostArray<Data,CPT>::operator = (DeviceArray<Data,CPT> & other)
  {
     Data* data = Array<Data>::data_;
     bool same =  (bool)data && data == other.cArray();
     if (same) {
        UTIL_CHECK(Array<Data>::capacity() == other.capacity());
     } else {
        FftwDRArray<Data>::associate(other); 
     }
     return *this;
  }

}
#endif
