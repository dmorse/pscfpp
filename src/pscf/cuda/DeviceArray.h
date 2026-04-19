#ifndef PSCF_DEVICE_ARRAY_H
#define PSCF_DEVICE_ARRAY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/ReferenceCounter.h>   // member
#include <util/misc/CountedReference.h>   // member

// Forward declarations
namespace Util {
   template <typename Data> class DArray;
}
namespace Pscf {
   class DeviceMemory;
}

namespace Pscf {

   using namespace Util;

   /**
   * Dynamic array on the GPU device with aligned data.
   *
   * This class wraps an aligned C array with elements of type Data that
   * is allocated in GPU device global memory.  All member functions may
   * be called from the CPU host, but this class does not offer access
   * to individual elements via the subscript operator, operator[].
   *
   * A DeviceArray can be in any of three states:
   *
   *   (1) Unallocated: In this state, there is no associated memory,
   *   so capacity() returns 0, while isAllocated(), isOwner() and
   *   isAssociated() all return false.
   *
   *   (2) A data owner: In this case, this object owns a block of memory
   *   that it is responsible for de-allocating. In this state, capacity()
   *   returns a positive integer value, isAllocated() and isOwner()
   *   return true, and isAssociated() returns false.
   *
   *   (3) A data user: In this case, this object has a pointer to a block
   *   of device memory that is owned by a different container. We describe
   *   this by saying this DeviceArray (the data user) is "associated"
   *   with or "references" memory owned by another container (the data
   *   owner). In this state, capacity() returns a positive
   *   value, isAllocated() and isAssociated() return true, and
   *   isOwner() returns false.
   *
   * Memory that is owned by a DeviceArray may be allocated and
   * deallocated by the allocate() and deallocate() member functions.
   * Memory owned by a DeviceArray is also dellocated when the
   * DeviceArray is destroyed.
   *
   * A DeviceArray<Data> object may be associated as a data user with:
   *
   *   - A slice of an array of Data objects owned by another
   *     DeviceArray<Data>, or
   *
   *   - A block of bare memory owned by a DeviceMemory container
   *
   * Associations with either type of data owner can be created and
   * destroyed by calling the associate and dissociate member functions
   * of the data user. Attempts to create an association cause an
   * Exception to be thrown if the requested slice of shared data
   * would exceed the bounds of the block owned by the data owner.
   *
   * A DeviceArray<Data> that serves as an owner of data that is used
   * by one or more other associated DeviceArray<Data> objects maintains
   * a count of how many associated objects refer to its data. This
   * counter is incremented when an association is created by a data
   * user and decremented when the user destroys the association.
   *
   * When a DeviceArray<Data> object is destroyed, any data that it
   * owns is automatically deallocated, and any association with
   * external data is destroyed.
   *
   * It is an error to deallocate memory owned by a DeviceArray or
   * DeviceMemory container that is still referred to by one or more
   * associated DeviceArray data users. Doing so is illegal because it
   * creates dangling pointers. Attempts to prematurely deallocate such
   * a shared memory block either cause an Exception to be thrown, if
   * it occurs within the deallocate function, or cause an error message
   * to be written to std::cout if the error occurs in the DeviceArray
   * destructor function. All references to such a shared data block
   * must thus be released, by calling the dissociate() member function
   * of each data user, before the data owner can be safely destroyed.
   *
   * \ingroup Pscf_Cuda_Containers_Module
   */
   template <typename Data>
   class DeviceArray
   {

   public:

      /**
      * Data type of each element.
      */
      typedef Data ValueType;

      /**
      * Default constructor.
      */
      DeviceArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      *
      * \param capacity  number of elements to allocate
      */
      DeviceArray(int capacity);

      /**
      * Copy constructor.
      *
      * \param other  DeviceArray<Data> to be copied (input)
      */
      DeviceArray(DeviceArray<Data> const & other);

      /**
      * Copy constructor, deep copy DArray<Data> from host to device.
      *
      * \param other  DArray<Data> to be copied (input)
      */
      DeviceArray(DArray<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DeviceArray();

      /**
      * Allocate the underlying C array on the device.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param capacity  number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the array is not allocated.
      * \throw Exception if this object does not own the array.
      */
      void deallocate();

      /**
      * Associate this object with a slice of a different DeviceArray.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param arr  parent array that owns the data
      * \param beginId  index in the parent array at which this array starts
      * \param capacity  number of elements associated with this container
      */
      void associate(DeviceArray<Data>& arr, int beginId, int capacity);

      /**
      * Associate this object with a DeviceMemory container.
      *
      * When a DeviceArray<Data> is associated with a DeviceMemory
      * container, the shared data block begins at the first byte of
      * the block owned by the DeviceMemory. The number of byes in such
      * a shared block is equal to capacity * sizeof(Data). This must
      * be less than or equal to the number of byes in the block owned
      * by the DeviceMemory, or this function will throw an Exception.
      *
      * \throw Exception if the array is already allocated
      * \throw Exception if the requested memory block is too large
      *
      * \param arr  DeviceMemory container that owns the data
      * \param capacity  number of elements of type Data
      */
      void associate(DeviceMemory& arr, int capacity);

      /**
      * Dissociate this object from an externally owned memory block.
      *
      * This function can destroy an association with either another
      * DeviceArray or with a DeviceMemory. After successful exit,
      * isAllocated() and isAssociated() will both return false.
      *
      * \throw Exception if this is not associated with external memory
      */
      void dissociate();

      /**
      * Associate a reference with reference counter for this container.
      *
      * This function should normally only be called within the associate
      * member function of a container that is referring to data owned
      * by this array, as part of the process of creating an association.
      * It should not be called by users in any other context.
      *
      * On exit, the ReferenceCounter owned by this DeviceMemory is
      * incremented by one and the CountedReference contains a pointer
      * to the ReferenceCounter. This pointer can later be used as part
      * of the process of destroying the association.
      *
      * \param reference  reference to be included by reference counter
      */
      void addReference(CountedReference& reference);

      /**
      * Assignment operator, assign from another DeviceArray<Data> array.
      *
      * Performs a deep copy, by copying values of all elements from
      * device memory to device memory.
      *
      * This function will allocate memory if this (LHS) array is not
      * allocated.  If this array is arleady allocated, it must have the
      * same capacity as the other (RHS) DeviceArray<Data>.
      *
      * \param other DeviceArray<Data> on rhs of assignent (input)
      */
      virtual
      DeviceArray<Data>& operator = (const DeviceArray<Data>& other);

      /**
      * Assignment operator, assignment from Util::DArray<Data> host array.
      *
      * Performs a deep copy from a RHS DArray<Data> host array to this 
      * LHS DeviceArray<Data>, by copying underlying C array from host 
      * memory to device memory.
      *
      * If this (LHS) DeviceArray<Data> is not allocated on entry, required
      * memory will be allocated before data is copied.  If this LHS object
      * is already allocated, it must have the same capacity as the RHS 
      * DArray<Data>.
      *
      * \param other  DArray<Data> on RHS of assignent (input)
      */
      virtual
      DeviceArray<Data>& operator = (const DArray<Data>& other);

      /**
      * Return array capacity.
      *
      * \return number of elements in array
      */
      int capacity() const;

      /**
      * Return true if the array has allocated data, false otherwise.
      *
      * A DeviceArray is considered allocated if it has non-null pointer
      * to an data block, which may either be a block that it owns or a
      * block owned by another container.
      */
      bool isAllocated() const;

      /**
      * Does this container own an allocated memory block?
      *
      * If isAllocated() is false, isOwner() must be false.
      */
      bool isOwner() const;

      /**
      * Is this container associated with a memory block it does not own?
      *
      * If isAllocated() is false, isAssociated() must be false.
      */
      bool isAssociated() const;

      /**
      * Return pointer to underlying C array.
      */
      Data* cArray();

      /**
      * Return const pointer to underlying C array.
      */
      Data const * cArray() const;

   protected:

      /// Pointer to a C array of Data elements on the GPU device.
      Data* dataPtr_;

      /// Allocated size (capacity) of the array.
      int capacity_;

      /// Counter for any arrays that use data owned by this.
      ReferenceCounter refCounter_;

      /// Reference to another array that owns memory used by this one.
      CountedReference ref_;

   };

   // Inline member function definitions

   /*
   * Return array capacity.
   */
   template <typename Data> inline
   int DeviceArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return true if this object has access to a memory block.
   */
   template <typename Data> inline
   bool DeviceArray<Data>::isAllocated() const
   {  return (bool) dataPtr_; }

   /*
   * Does this object own data?
   */
   template <typename Data> inline
   bool DeviceArray<Data>::isOwner() const
   {  return ((bool) dataPtr_ && !ref_.isAssociated()); }

   /*
   * Is this object associated with data it does not own?
   */
   template <typename Data> inline
   bool DeviceArray<Data>::isAssociated() const
   {  return ((bool) dataPtr_ && ref_.isAssociated()); }

}

#include "DeviceMemory.h"
#include "cudaErrorCheck.h"
#include <util/containers/DArray.h>
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Default constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray()
    : dataPtr_(nullptr),
      capacity_(0)
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(int capacity)
    : dataPtr_(nullptr),
      capacity_(0)
   {  allocate(capacity); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(const DeviceArray<Data>& other)
    : dataPtr_(nullptr),
      capacity_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other array must be allocated.");
      }

      allocate(other.capacity_);
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyDeviceToDevice) );
   }

   /*
   * Copy constructor, deep copy from host DArray<Data>
   *
   * Allocates new memory and copies all elements by value.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(const DArray<Data>& other)
    : dataPtr_(nullptr),
      capacity_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other array must be allocated.");
      }

      allocate(other.capacity());
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyHostToDevice) );
   }

   /*
   * Destructor.
   */
   template <typename Data>
   DeviceArray<Data>::~DeviceArray()
   {
      if (isOwner()) {
         if (refCounter_.hasRefs()) {
            std::cout
                << "Error: Destroying DeviceArray that is referred"
                << "to by at least one CountedReference"
                << std::endl;
         }
         cudaFree(dataPtr_);
      }
      if (ref_.isAssociated()) {
         ref_.dissociate();
      }
   }

   /*
   * Allocate device memory, with this object then owns.
   */
   template <typename Data>
   void DeviceArray<Data>::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate an array");
      }
      if (ref_.isAssociated()) {
         UTIL_THROW("Attempt to allocate array that uses external data");
      }

      cudaErrorCheck(cudaMalloc((void**) &dataPtr_,
                     capacity * sizeof(Data)));
      capacity_ = capacity;
   }

   /*
   * Deallocate device memory owned by this object, if any.
   */
   template <typename Data>
   void DeviceArray<Data>::deallocate()
   {
      UTIL_CHECK(dataPtr_);
      UTIL_CHECK(!ref_.isAssociated());
      UTIL_CHECK(!refCounter_.hasRefs());
      cudaErrorCheck( cudaFree(dataPtr_) );
      dataPtr_ = nullptr;
      capacity_ = 0;
   }

   /*
   * Associate this object with memory owned by a different DeviceArray.
   */
   template <typename Data>
   void DeviceArray<Data>::associate(DeviceArray<Data>& arr,
                                     int beginId, int capacity)
   {
      UTIL_CHECK(arr.isAllocated());
      UTIL_CHECK(arr.isOwner());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(beginId + capacity <= arr.capacity());
      UTIL_CHECK(!dataPtr_);
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and capacity
      dataPtr_ = arr.cArray() + beginId;
      capacity_ = capacity;

      // Associate ref_ member with reference counter of data owner
      arr.addReference(ref_);
   }

   /*
   * Associate this object with memory owned by a DeviceMemory object.
   */
   template <typename Data>
   void DeviceArray<Data>::associate(DeviceMemory& arr, int capacity)
   {
      UTIL_CHECK(arr.isAllocated());
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(capacity * sizeof(Data) <= arr.capacity());
      UTIL_CHECK(!dataPtr_);
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and capacity
      dataPtr_ = (Data *) arr.cArray();
      capacity_ = capacity;

      // Associate ref_ member with reference counter of data owner
      arr.addReference(ref_);
   }

   /*
   * Dissociate this object from external device memory
   */
   template <typename Data>
   void DeviceArray<Data>::dissociate()
   {
      UTIL_CHECK(dataPtr_);
      UTIL_CHECK(ref_.isAssociated());

      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
      ref_.dissociate();
   }

   /*
   * Associate an external CountedReference with the ReferenceCounter.
   *
   * This function should called to create an association in which
   * this object is the data owner. The CountedReference parameter
   * should be owned by the data user.
   */
   template <typename Data>
   void DeviceArray<Data>::addReference(CountedReference& ref)
   {  ref.associate(refCounter_); }

   /*
   * Assignment from another DeviceArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>&
   DeviceArray<Data>::operator = (const DeviceArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - RHS DeviceArray<Data> must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("Other DeviceArray<Data> must be allocated.");
      }

      // If this is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Require equal capacity values
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy elements
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyDeviceToDevice) );

      return *this;
   }

   /*
   * Assignment from RHS Util::DArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>&
   DeviceArray<Data>::operator = (const Util::DArray<Data>& other)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DArray<Data> must be allocated.");
      }

      // Allocate this if necessary
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Require equal capacity values
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy elements
      cudaErrorCheck( cudaMemcpy(dataPtr_, other.cArray(),
                                 capacity_ * sizeof(Data),
                                 cudaMemcpyHostToDevice) );

      return *this;
   }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   Data* DeviceArray<Data>::cArray()
   {
      UTIL_CHECK(dataPtr_);
      return dataPtr_;
   }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   const Data* DeviceArray<Data>::cArray() const
   {
      UTIL_CHECK(dataPtr_);
      return dataPtr_;
   }

} // namespace Pscf
#endif
