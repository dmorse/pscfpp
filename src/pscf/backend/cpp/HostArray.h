#ifndef PSCF_HOST_ARRAY_CP_H
#define PSCF_HOST_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>       // base class
#include <util/misc/CountedReference.h>  // member

// Forward declarations
namespace Pscf {
   template <typename Data, typename T> class DeviceArray;
   class CPT;
}

namespace Pscf {

   using namespace Util;

   // Declare primary template
   template <typename Data, typename T> class HostArray;

   /**
   * Psuedo-"host" array that provides read-write access.
   *
   * This class template partial specialization may be used in template
   * code in which data is assigned to or from a psuedo-device array to 
   * or from a pseudo-host array, in cases when read-write access to the 
   * host array is required. A ConstHostArray<Data> should be used for
   * device-to-host data transfer in which only read access is required
   * on the host, and/or the device array is declared const.
   *
   * Construction or assignment (operator =)from a DeviceArray<Data,CPT>
   * and the "associate" member function all create a shallow read-only 
   * copy of a C array that is owned by an associated device array. The
   * assignment operator does nothing if the desired association already
   * exists. Each of these functions takes a non-const reference to the
   * DeviceArray<Data,CPT> as a parameter, which is not a const reference
   * because they all provide write access to the underlying shared array.
   *
   * The lifetime of an association with a device array must not be 
   * allowed to extend beyond the function in which the association is
   * created. This is necessary to guarantee that no associations remain
   * when the device array is deallocated or destroyed.  Such an 
   * association can be released by calling the dissociate function of 
   * the host array, or will be released upon destruction of the host 
   * array.  A reference counting system is used to detect and report 
   * erroneous de-allocation of a device array when it is still referred 
   * to by one or more other host arrays.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class HostArray<Data, CPT> : public Array<Data>
   {

   public:

      /**
      * Default constructor.
      */
      HostArray();

      // Prohibit copy construction from another HostArray.
      HostArray(HostArray<Data,CPT> const & other) = delete;

      /**
      * Conversion constructor from a device array.
      *
      * Create an association (shallow copy) with memory owned by the
      * pre-existing device array, thus creating a shallow copy that
      * points to the same memory. This is equivalent to default 
      * construction followed by the associate function.
      *
      * \throw Exception if the other device array is not allocated
      *
      * \param other  associated device array
      */
      HostArray(DeviceArray<Data,CPT> & other);

      /**
      * Destructor.
      */
      ~HostArray();

      // Prohibit assignment from another HostArray.
      HostArray<Data,CPT>&
      operator = (HostArray<Data,CPT> const&) = delete;

      /**
      * Associate this object with a device array.
      *
      * Associates this object with a device array, thus creating a
      * shallow copy that points to the same memory. After successful 
      * return, isAllocated() will return true.
      *
      * \throw Exception if source array is not allocated on entry
      * \throw Exception if this array is already associated
      *
      * \param source  array that owns the data
      */
      void associate(DeviceArray<Data,CPT> & other);

      /**
      * Pseudo-assignment from a Device array.
      *
      * If this is already associated with the other device array that
      * is passed as a argument, this function does nothing and returns.
      * Otherwise, if this object is not associated with any data source,
      * the assignment operator creates an association with the specified
      * device array, thus creating a shallow copy.
      *
      * \throw Exception if RHS device array is not allocated
      * \throw Exception if this is associated with any other data source
      *
      * \param other  device array on RHS of assignment 
      */
      HostArray<Data,CPT>& 
      operator = (DeviceArray<Data,CPT> & other);

      /**
      * Release association with a device array.
      *
      * Upon successful return, isAllocated() will return false.
      *
      * \throw Exception if this is not associated with another array.
      */
      void dissociate();

   private:

      /// Reference to a device array that has an array referenced by this.
      CountedReference ref_;

      // Hide inherited protected members by making them private
      using Array<Data>::data_;
      using Array<Data>::capacity_;

   };

} // namespace Pscf

#include <pscf/backend/cpp/DeviceArray.h>
#include <util/global.h>

namespace Pscf {

   // Member functions

   /*
   * Default constructor.
   */
   template <typename Data>
   HostArray<Data,CPT>::HostArray()
    : Array<Data>(),
      ref_()
   {}

   /*
   * Conversion construction from a device array.
   *
   * Creates an association with a DeviceArray, which must be allocated.
   */
   template <typename Data>
   HostArray<Data,CPT>::HostArray(DeviceArray<Data,CPT> & other)
    : Array<Data>(),
      ref_()
   {  associate(other); }

   /*
   * Destructor.
   */
   template <typename Data>
   HostArray<Data,CPT>::~HostArray()
   {
      if (ref_.isAssociated()) {
         ref_.dissociate(); // decrements counter of source array
      }
   }

   /*
   * Associate this object with a device array.
   */
   template <typename Data>
   void 
   HostArray<Data,CPT>::associate(DeviceArray<Data,CPT> & source)
   {
      UTIL_CHECK(source.isAllocated());
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and size
      data_ = source.cArray();
      capacity_ = source.capacity();

      // Note: The const_cast to non-const pointer is permissible because
      // the public interface of the Array<Data> base class is designed
      // prevent modification of the values of array elements.

      // Associate ReferenceCounter of the source array with this 
      ref_.associate(source);

      // On exit, the ReferenceCounter sub-object of the data source is 
      // incremented and the ref_ CountedReference member variable of 
      // this object holds a pointer to that ReferenceCounter.
   }

   /*
   * Assignment from a DeviceArray.
   *
   * Creates an association with a DeviceArray, or do nothing if the
   * association already exists.
   */
   template <typename Data>
   HostArray<Data,CPT>&
   HostArray<Data,CPT>::operator = (DeviceArray<Data,CPT>& other)
   {
      UTIL_CHECK(other.isAllocated());
      if (!ref_.isAssociated()) {
         associate(other);
      } else {
         Data * data = Array<Data>::cArray();
         UTIL_CHECK((bool)data);
         UTIL_CHECK(other.cArray() == data);
         UTIL_CHECK(other.capacity() == Array<Data>::capacity());
      }
      return *this;
   }

   /*
   * Release pre-existing association with a device array.
   */
   template <typename Data>
   void HostArray<Data,CPT>::dissociate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(ref_.isAssociated());

      data_ = nullptr;
      capacity_ = 0;
      ref_.dissociate(); // decrements reference counter of source array
   }

}
#endif
