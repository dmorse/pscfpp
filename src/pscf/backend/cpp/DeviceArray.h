#ifndef PSCF_DEVICE_ARRAY_CP_H
#define PSCF_DEVICE_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FftwDRArray.h"         // base class
#include "CPT.h"                 // template argument

// Forward declaration
namespace Pscf {
   template <typename Data, typename T> class HostArray;
}

namespace Pscf {

   using namespace Util;

   // Declare primary template
   template <typename Data, typename T> class DeviceArray;

   /**
   * Pseudo "device" array for use with C++ backend.
   *
   * Derived from FftwDRArray, and largely equivalent. 
   * 
   * The main difference from the base class is the addition of assignment 
   * (operator = ) from a HostArray<Data,CTP> to a DeviceArray<Data,CTP>.
   * This operator checks for the existence of an association between the
   * device array (LHS) and host array (RHS) in which both point to the
   * same memory, throws an Exception if such an association does not
   * already exist, or does nothing if it does. Reason for this usage is
   * discussed below.
   *
   * Host to device data transfer:
   * 
   * Usage for transfer of data from an instance of HostArray<Data,CPT>
   * to an instance of DeviceArray<Data,CPT> is discussed in the class
   * documentation for the HostArray class template. The most explicit
   * version of this pattern involves the following steps:
   *
   *    - Construct and allocate a device array (denoted by dArray)
   *
   *    - Default construct a host array (denoted by hArray)
   * 
   *    - Create an association: hArray.associate(dArray)
   *
   *    - Initialize elements of hArray
   *
   *    - Assign from host to device: dArray = hArray;
   *
   *    - Destroy the association: hArray.dissociate()
   *
   *    - Deallocate or destroy the device array dArray
   *
   * Comments:
   *
   *    - The operator v = u that assigns from host u to device v merely
   *      checks that an association already exists and does nothing if
   *      it does. This is because the association must exist before the 
   *      data is initialized on the host array, which must occur before
   *      the assignment operator is invoked. In analogous GPU code, the
   *      assignment operator would instead actually transfer data from 
   *      CPU to GPU memory.
   *
   *    - The dissociate member function of the host array class destroys
   *      an association in CPU code (T=CPT) but does nothing in analogous 
   *      GPU code (T=CUT).
   *
   *    - The host array must never be used to modify data after the 
   *      dissociate function is invoked. To enforce this, dissociate
   *      should usually be called immediately after the assignment 
   *      operator. 
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class DeviceArray<Data, CPT> : public FftwDRArray<Data>
   {

   public:

      using typename FftwDRArray<Data>::ValueType;

      /**
      * Backend identifier class typename aliase.
      */
      using BackendIdClass = CPT;

      // Default constructor
      DeviceArray() = default;

      /**
      * Allocating constructor.
      *
      * \param capacity number of elements to allocate
      */
      DeviceArray(int capacity);

      // Copy constructor.
      DeviceArray(DeviceArray<Data,CPT> const & other) = default;

      // Destructor.
      ~DeviceArray() override = default;

      // Assignment.
      DeviceArray<Data,CPT>& 
      operator = (DeviceArray<Data,CPT> const&) = default;

      /**
      * Pseudo-assignment from a host array.
      *
      * This function performs a sanity check and does nothing if the 
      * check is passed. It checks that both LHS and RHS are allocated
      * and refer to the to the same underlying C array on entry, and
      * throws an Exception if this is not the case. 
      *
      * Rationale: Since the host array acts as a shallow copy of this 
      * device array, the association must have been created before data 
      * was initialized on the host array, which must occurs before this
      * transaction is finalized by the assignment operator.
      *
      * \throw Exception if this is not allocated
      * \throw Exception if other array is not allocated
      * \throw Exception if this and other do not point to the same memory
      *
      * \param other  array container on RHS of assigment (input)
      */
      DeviceArray<Data,CPT>& operator = (HostArray<Data,CPT> const & other);

      // Inherited public function (to prevent hiding)
      using FftwDRArray<Data>::operator =;

   };

   // Explicit instantiation declarations
   extern template class DeviceArray<double,CPT>; 
   extern template class DeviceArray<fftw_complex,CPT>; 

} // namespace Pscf

#include "HostArray.h"
namespace Pscf {

   /*
   * Allocating constructor.
   */
   template <typename Data>
   DeviceArray<Data,CPT>::DeviceArray(int capacity)
    : FftwDRArray<Data>(capacity)
   {}

   /*
   * Check for a pre-existing association with a HostArray.
   */
   template <typename Data>
   DeviceArray<Data,CPT>& 
   DeviceArray<Data,CPT>::operator = (HostArray<Data,CPT> const & other)
   {
      UTIL_CHECK(Array<Data>::isAllocated());
      UTIL_CHECK(other.isAllocated());
      UTIL_CHECK(Array<Data>::cArray() == other.cArray());
      UTIL_CHECK(Array<Data>::capacity() == other.capacity());
      return *this;
   }

}
#endif
