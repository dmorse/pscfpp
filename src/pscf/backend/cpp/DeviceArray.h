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
   * In order to transfer data from HostArray<Data,CPT> u to a
   * DeviceArray<Data,CPT> v one must:
   *
   *    - Allocate the device array
   *
   *    - Create an association either by using the conversion constructor
   *      HostArray<Data,CPT> u(v) to create the host array or by calling 
   *      the host array associate function u.associate(v).
   *
   *    - Initialize data on the host array u
   *
   *    - Invoke the assignment operator: v = u;
   *
   *    - Destroy the association, either by explicitly invoking the
   *      dissociate function u.dissociate() on the host array or by 
   *      allowing a host array that is a local object to be destroyed
   *      when it goes out of scope
   *
   * Comments:
   *
   *    - The operator v = u that assigns from host u to device v merely
   *      checks that an association already exists and does nothing if
   *      it does. This is because the association must exist before the 
   *      data is initialized on the host array, which must occur before
   *      the assignment operator is invoked. In analogous GPU code, the
   *      assignment operator would instead transfer data from CPU to GPU
   *      memory.
   *
   *    - The lifetime of association between host and device arrays 
   *      may never be allowed to extend beyond the function in which it 
   *      was created. This is necessary to guarantee that an association 
   *      will never still exist when the shared array is deallocated.
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
      ~DeviceArray() = default;

      // Assignment.
      DeviceArray<Data,CPT>& 
      operator = (DeviceArray<Data,CPT> const&) = default;

      /**
      * Pseudo-assignment from a host array.
      *
      * This functions checks that both arrays are allocated and refer to
      * the same underlying C array, and throws an Exception if this is 
      * not the case. Rationale: Since the host array acts as a shallow
      * copy of this device array, the association must have been created
      * before data was initialized on the host array, which occurs before
      * the transaction is finalized by the assignment operator.
      *
      * This function also dissociates the host array from this device array. 
      * Rationale: In assignment from host to device, the assignment operator 
      * marks the end of the transaction. In the corresponding code with T=CUT,
      * the persistent data on the finalized is finalized by this operator. 
      * In this specialization, with T=CPT, leaving an association live would 
      * allow later modification of data on the host to change data in the 
      * device array, which is inconsistent with the behavior of the CUDA
      * specialization. 
      *
      * \throw Exception if this is not allocated
      * \throw Exception if other array is not allocated
      * \throw Exception if this and other do not point to the same memory
      *
      * \param other  array container on RHS of assigment (input)
      */
      DeviceArray<Data,CPT>& operator = (HostArray<Data,CPT> & other);

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
   DeviceArray<Data,CPT>::operator = (HostArray<Data,CPT> & other)
   {
      UTIL_CHECK(Array<Data>::isAllocated());
      UTIL_CHECK(other.isAllocated());
      UTIL_CHECK(Array<Data>::cArray() == other.cArray());
      UTIL_CHECK(Array<Data>::capacity() == other.capacity());
      other.dissociate();
      return *this;
   }

}
#endif
