#ifndef PSCF_HOST_RECV_ARRAY_CU_H
#define PSCF_HOST_RECV_ARRAY_CU_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>       // template argument
#include <util/containers/ConstArray.h>  // base class template

// Forward declarations
namespace Util {
   template <typename Data> class Array;
}
namespace Pscf {
   template <typename Data, typename T> class HostRecvArray;
   template <typename Data, typename T> class DeviceArray;
}

namespace Pscf {

   using namespace Util;

   /**
   * Host array to receive host-to-device data transfers.
   *
   * This class template should be used in backend-independent template
   * code in which data is copied from a GPU device array to a CPU host 
   * array.  This template is derived from ConstArray<Data>, and thus 
   * provides read-only access to the data.
   *
   * Conversion construction or assignment (operator =) from an instance 
   * of DeviceArray<Data,CUT> to a HostRecvArray<Data,CUT> allocates 
   * the array if it is not allocated and then copies all elements of 
   * the array from global GPU device memory to host memory. Conversion
   * construction from a device array is equivalent to default 
   * construction followed by assignment.
   *
   * <b> Usage: </b>
   *
   * Typical usage is shown below for backend-indepent template code 
   * for device-to-host transfer from a longer lived instance of 
   * DeviceArray<Data,CUT> named deviceArray to a local object that is
   * an instance of HostRecvArray<Data,CUT> named hostArray. The type
   * type of each array element is denoted by Data.
   * \code
   *    HostRecvArray<Data,CUT> hostArray;
   *    hostArray = deviceArray
   *
   *    // Read data from the hostArray 
   *
   *    hostArray.dissociate();
   * \endcode
   * If the hostArray is a local object that is defined with the 
   * function in which it is used, the default constructor and 
   * assignment operations may instead be combined into a single 
   * call of the conversion constructor, giving the shorter version
   * \code
   *    HostRecvArray<Data,CUT> hostArray(deviceArray)
   *
   *    // Read data from the hostArray 
   *
   *    hostArray.dissociate();
   * \endcode
   * In this case, the call of the dissociate function may also be
   * omitted, since an equivalent operation will be performed by
   * the destructor when the function ends and the hostArray goes 
   * out of scope. 
   *
   * Comments:
   * 
   *  - In this specialization for a CUDA backend (T=CUT), the assignment
   *    operator and conversion constructor each perform a deep copy of 
   *    data from GPU device memory to CPU host memory. 
   *
   *  - This class template provides read-only access to data that is
   *    copied from the GPU device. If the copied data must be modified,
   *    it must first be copied to a CPU container that allows this.
   *
   *  - In this specialization for a CUDA backend (T=CUT), the dissociate 
   *    function does nothing. It is used in backend-independent code to
   *    maintain compatibility with the syntax used by the specialization
   *    for a C++ backend (T=CPT), for which the dissociate function
   *    releases an association with the device array.
   *
   *  - Any call to the dissociate function must occur after all statements 
   *    that read data from the host array. 
   *
   * \ingroup Pscf_Backend_Cuda_Module
   */
   template <typename Data>
   class HostRecvArray<Data,CUT> : public ConstArray<Data>
   {

   public:

      // Public types

      /**
      * Alias for data type of each array element.
      */
      using ValueType = Data;

      /**
      * Alias for backend identifier class.
      */
      using BackendIdClass = CUT;

      // Public member functions

      /**
      * Default constructor.
      */
      HostRecvArray() = default;

      /**
      * Allocating constructor.
      *
      * \param capacity  number of elements to allocate
      */
      HostRecvArray(int capacity);

      // Default copy constructor (delete).
      HostRecvArray(HostRecvArray<Data,CUT> const &) = delete;

      /**
      * Conversion constructor (copies from device to host).
      *
      * Conversion construction is equivalent to default construction 
      * followed by assignment. This function first allocates this array
      * with a capacity equal to that of the other device array, and then
      * performs a deep copy of the underlying array from GPU device 
      * memory to CPU host memory.
      *
      * \throw Exception if the RHS array is not allocated on entry
      *
      * \param other DeviceArray<Data,CUT> to be copied (input)
      */
      HostRecvArray(DeviceArray<Data,CUT> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~HostRecvArray();

      // Default assignment (delete)
      HostRecvArray<Data,CUT>& 
      operator = (HostRecvArray<Data,CUT> const &) = delete;

      /**
      * Assignment from an Array<Data> container.
      *
      * Performs a deep copy, by copying values of all elements of an
      * Array<Data> container. If this (LHS) array is already allocated
      * on entry, it must have the same capacity as the other (RHS) array.
      * If this LHS array is not allocated on entry, required memory is
      * allocated before copying values.
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if arrays are allocated with unequal capacities
      *
      * \param other  array container on RHS of assigment (input)
      */
      HostRecvArray<Data,CUT>& operator = (Array<Data> const & other);

      /**
      * Assignment from a DeviceArray<Data,CUT>.
      *
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this 
      * LHS HostRecvArray<D>, by first allocating if necessary, and then
      * copying the underlying data from device memory to host memory.
      *
      * \throw Exception if the RHS array is not allocated on entry
      * \throw Exception if LHS and RHS have unequal nonzero capacities
      *
      * \param other  DeviceArray<Data,CUT> on RHS of assignment (input)
      */
      HostRecvArray<Data,CUT>&
      operator = (DeviceArray<Data,CUT> const & other);

      /**
      * Copy a slice of data from a larger DeviceArray into this array.
      *
      * This function populate this HostRecvArray with data from a slice
      * of a DeviceArray. The size of the slice is equal to the capacity of
      * this LHS HostRecvArray, so that this entire array will be fully
      * populated. The position of the beginning of the slice within the 
      * RHS DeviceArray is indicated by input parameter beginId. The
      * capacity of the RHS DeviceArray must thus be greater than or equal
      * to the sum of beginId and the capacity of this LHS HostRecvArray.
      *
      * \throw Exception if RHS array is not allocated
      * \throw Exception if LHS array is not allocated
      * \throw Exception if beginId is negative
      * \throw Exception if the slice extends beyond the end of the RHS
      *
      * \param other DeviceArray<Data,CUT>  object from which to copy slice
      * \param beginId  index of other array at which slice begins
      */
      void copySlice(DeviceArray<Data,CUT> const & other, int beginId);

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the HostRecvArray is already allocated
      *
      * \param capacity  number of elements to allocate
      */
      void allocate(int capacity);

      /**
      * Deallocate the underlying C array.
      *
      * \throw Exception if the HostRecvArray is not allocated
      */
      void deallocate();

      #if 0
      /**
      * Reallocate the underlying C array if necessary.
      *
      * The array is reallocated and copied to a new location if the new
      * capacity, given by the capacity parameter, is greater than the
      * existing array capacity. Nothing is done if the new and old
      * capacities are equal. An Exception is thrown if the new capacity
      * is less than the old capacity.
      *
      * \param capacity  number of elements for which to allocate space
      */
      void reallocate(int capacity);
      #endif

      /**
      * Release this array after finishing all reading.
      *
      * In template code that must work with either CPU or GPU, this should
      * be called after all code that requires read access to this array.
      * In this GPU specialization (T=CUT), this function does nothing. 
      */
      void dissociate()
      {}

      /**
      * Serialize a HostRecvArray to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      // Inherited public member functions (selected)
      using ConstArray<Data>::operator [];
      using ConstArray<Data>::cArray;
      using ConstArray<Data>::isAllocated;
      using ConstArray<Data>::capacity;
      using ConstArray<Data>::size;

   private:

      // Hide inherited protected member variables
      using ConstArray<Data>::data_;
      using ConstArray<Data>::capacity_;

   };

   // Explicit instantiation declarations
   extern template class HostRecvArray<cudaReal,CUT>;
   extern template class HostRecvArray<cudaComplex,CUT>;

} // namespace Pscf

#include "DeviceArray.h"
#include "cudaErrorCheck.h"
#include <util/containers/Array.h>
#include <util/misc/Memory.h>
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Allocating constructor.
   */
   template <typename Data>
   HostRecvArray<Data,CUT>::HostRecvArray(int capacity)
    : ConstArray<Data>()
   {  allocate(capacity); }

   /*
   * Conversion constructor from DeviceArray<Data,CUT>.
   *
   * Allocate and perform a deep copy from device to host.
   */
   template <typename Data>
   HostRecvArray<Data,CUT>::HostRecvArray(
      DeviceArray<Data,CUT> const& other)
    : ConstArray<Data>()
   {  (*this) = other; }

   /*
   * Destructor.
   */
   template <typename Data>
   HostRecvArray<Data,CUT>::~HostRecvArray()
   {
      if (isAllocated()) {
         try {
            Memory::deallocate<Data>(data_, capacity_);
         } catch (...) {
            data_ = nullptr;
            std::cout << "Exception in HostRecvArray destructor";
         }
         capacity_ = 0;
      }
   }

   /*
   * Assignment from an Array<Data> (deep copy on host).
   */
   template <typename Data>
   HostRecvArray<Data,CUT>& 
   HostRecvArray<Data,CUT>::operator = (Array<Data> const & other)
   {
      // Check for pseudo self assignment
      if (dynamic_cast< Array<Data>* >(this) == &other) return *this;

      // Preconditions - require that other (RHS) array is allocated
      UTIL_CHECK(other.isAllocated());
      UTIL_CHECK(other.capacity() > 0);

      // If this LHS array is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Require equal capacities
      UTIL_CHECK (capacity_ == other.capacity());

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Assignment from a DeviceArray<Data,CUT> (device-to-host copy).
   *
   * Allocate if necessary, and perform a deep copy from device to host.
   */
   template <typename Data>
   HostRecvArray<Data,CUT>&
   HostRecvArray<Data,CUT>::operator = (
                                  DeviceArray<Data,CUT> const & other)
   {
      // Preconditions on other array
      UTIL_CHECK(other.isAllocated());
      const int n = other.capacity();
      UTIL_CHECK(n > 0);

      // Allocate this if not allocated previously. 
      if (!isAllocated()) {
         allocate(n);
      }
      UTIL_CHECK(capacity() == n);
      // If this was allocated with capacity() == n, do nothing.
      // If this was allocated with capacity() != n, throw Exception.

      // Copy all elements
      cudaErrorCheck(
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray(),
                    capacity() * sizeof(Data),
                    cudaMemcpyDeviceToHost)
      );

      return *this;
   }

   /*
   * Copy a slice of the data from a larger DeviceArray into this array.
   */
   template <typename Data>
   void HostRecvArray<Data,CUT>::copySlice(
                                    DeviceArray<Data,CUT> const & other,
                                    int beginId)
   {
      // Preconditions
      UTIL_CHECK (other.isAllocated());
      UTIL_CHECK(isAllocated());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity() + beginId <= other.capacity());

      // Copy all elements
      cudaErrorCheck(
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray() + beginId,
                    capacity() * sizeof(Data),
                    cudaMemcpyDeviceToHost)
      );

   }

   /*
   * Allocate the underlying C array.
   */
   template <typename Data>
   void HostRecvArray<Data,CUT>::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a HostRecvArray");
      }
      Memory::allocate<Data>(data_, capacity);
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   */
   template <typename Data>
   void HostRecvArray<Data,CUT>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      Memory::deallocate<Data>(data_, capacity_);
      capacity_ = 0;
   }

   #if 0
   /*
   * Reallocate the underlying C array, if necessary.
   */
   template <typename Data>
   void HostRecvArray<Data,CUT>::reallocate(int capacity)
   {
      UTIL_CHECK(capacity >= 0);
      if (capacity == capacity_) return;

      UTIL_CHECK(capacity > capacity_);
      if (isAllocated()) {
         Memory::reallocate<Data>(data_, capacity_, capacity);
      } else {
         Memory::allocate<Data>(data_, capacity);
      }
      capacity_ = capacity;
   }
   #endif

   /*
   * Serialize a HostRecvArray to/from an Archive.
   */
   template <typename Data>
   template <class Archive>
   void HostRecvArray<Data,CUT>::serialize(Archive& ar, 
                                            const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent HostRecvArray capacities");
            }
         }
      }
      if (isAllocated()) {
         for (int i = 0; i < capacity_; ++i) {
            ar & data_[i];
         }
      }
   }

}
#endif
