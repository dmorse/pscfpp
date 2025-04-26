#ifndef PSCF_DEVICE_ARRAY_H
#define PSCF_DEVICE_ARRAY_H

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "cudaErrorCheck.h"
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   using namespace Util;

   // Forward declaration of analogous container for data on host.
   template <typename Data> class HostDArray;

   /**
   * Dynamic array on the GPU device with aligned data.
   *
   * This class wraps an aligned C array with elements of type Data that 
   * is allocated in GPU device global memory.  All member functions may 
   * be called from the CPU host, but this class does not offer access 
   * to individual elements via the subscript operator, operator[].
   * 
   * A DeviceArray can have one of two different relationships with its
   * underlying data. In case 1, the data are owned by this object, so
   * the allocation and destruction of the C array are performed by this
   * object. In case 2, this object is instead "associated" with a slice
   * of an array owned by a different DeviceArray. If this is the case, 
   * this object is not responsible for allocation or destruction of the 
   * underlying C array, and merely acts as a reference to the slice of
   * the other array to which it is associated. 
   *
   * \ingroup Pscf_Cuda_Module
   */
   template <typename Data>
   class DeviceArray
   {

   public:

      /**
      * Data type of each element.
      */
      typedef Data ElementType;

      /**
      * Default constructor.
      */
      DeviceArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity number of elements to allocate 
      */
      DeviceArray(int capacity);

      /**
      * Copy constructor.
      * 
      * \param other DeviceArray<Data> to be copied (input)
      */
      DeviceArray(DeviceArray<Data> const & other);

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
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Associate this object with a slice of a different DeviceArray.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param arr parent array that owns the data
      * \param beginId index in the parent array at which this array starts
      * \param capacity number of elements to allocate 
      */
      void associate(DeviceArray<Data>& arr, int beginId, int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the array is not allocated.
      * \throw Exception if this object does not own the array.
      */
      void deallocate();

      /**
      * Dissociate this object from the associated array.
      *
      * \throw Exception if this object is not associated with an array.
      */
      void dissociate();

      /**
      * Assignment operator, assign from another DeviceArray<Data> array.
      *  
      * Performs a deep copy, by copying values of all elements from 
      * device memory to device memory.
      *
      * This function will allocate memory if this (LHS) array is not 
      * allocated.  If this is allocated, it must have the same 
      * dimensions as the RHS DeviceArray<Data>.
      *
      * \param other DeviceArray<Data> on rhs of assignent (input)
      */
      virtual 
      DeviceArray<Data>& operator = (const DeviceArray<Data>& other);

      /**
      * Assignment operator, assignment from HostDArray<Data> host array.
      *
      * Performs a deep copy from a RHS HostDArray<Data> host array to 
      * this LHS DeviceArray<Data> device array, by copying underlying 
      * C array from host memory to device memory.
      *
      * This function will allocate memory if this (LHS) 
      * DeviceArray<Data> is not allocated.  If this is allocated, it 
      * must have the same dimensions as the RHS HostDArray<Data>.
      *
      * \param other HostDArray<Data> on RHS of assignent (input)
      */
      virtual 
      DeviceArray<Data>& operator = (const HostDArray<Data>& other);

      /**
      * Return allocated capacity.
      *
      * \return Number of elements allocated in array
      */
      int capacity() const;

      /**
      * Return true if the array has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Does this object own the underlying data? 
      */
      bool isOwner() const;

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

      /**
      * Pointer to array that owns this data, if isOwner_ == false.
      * 
      * Used to check that parent array has not been deallocated and/or
      * reallocated.
      */
      DeviceArray<Data> const * ownerPtr_;

      /**
      * Capacity of parent array, if isOwner_ == false.
      * 
      * Used to check that parent array has not been deallocated and/or
      * reallocated.
      */
      int ownerCapacity_;

      /**
      * Const pointer to parent array, if isOwner_ == false
      * 
      * Used to check that parent array has not been deallocated and/or
      * reallocated.
      */
      Data const * ownerDataPtr_;
   };

   /*
   * Default constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray()
    : dataPtr_(nullptr),
      capacity_(0),
      ownerPtr_(nullptr),
      ownerCapacity_(0),
      ownerDataPtr_(nullptr)
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(int capacity)
    : dataPtr_(nullptr),
      capacity_(0),
      ownerPtr_(nullptr),
      ownerCapacity_(0),
      ownerDataPtr_(nullptr)
   {  allocate(capacity); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <typename Data>
   DeviceArray<Data>::DeviceArray(const DeviceArray<Data>& other)
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
   * Destructor.
   */
   template <typename Data>
   DeviceArray<Data>::~DeviceArray()
   {
      if (isAllocated() && isOwner()) {
         cudaErrorCheck( cudaFree(dataPtr_) );
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   template <typename Data>
   void DeviceArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate an array");
      }
      if (ownerPtr_) {
         UTIL_THROW("Attempt to allocate array already associated with data");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      cudaErrorCheck( cudaMalloc((void**) &dataPtr_, capacity * sizeof(Data)) );
      capacity_ = capacity;

      ownerPtr_ = nullptr;
   }

   /*
   * Associate this object with a slice of a different DeviceArray.
   */
   template <typename Data>
   void DeviceArray<Data>::associate(DeviceArray<Data>& arr, int beginId, 
                                     int capacity)
   {
      UTIL_CHECK(arr.isAllocated());
      UTIL_CHECK(arr.isOwner());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(beginId + capacity <= arr.capacity());
      
      if (isAllocated()) {
         UTIL_THROW("Attempt to associate an already-allocated array.");
      }
      
      dataPtr_ = arr.cArray() + beginId;
      capacity_ = capacity;
      ownerPtr_ = &arr;
      ownerCapacity_ = arr.capacity();
      ownerDataPtr_ = arr.cArray();
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this array is not allocated.
   */
   template <typename Data>
   void DeviceArray<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      if (!isOwner()) {
         UTIL_THROW("Cannot deallocate, data not owned by this object.");
      }
      cudaErrorCheck( cudaFree(dataPtr_) );
      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
   }

   /*
   * Dissociate this object from the associated array.
   *
   * Throw Exception if this object is not associated with another array.
   */
   template <typename Data>
   void DeviceArray<Data>::dissociate()
   {
      if (isOwner()) {
         UTIL_THROW("Cannot dissociate: this object owns the array.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Cannot dissociate: no associated data found.");
      }
      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
      ownerPtr_ = nullptr; // reset to null
      ownerCapacity_ = 0;
      ownerDataPtr_ = nullptr;
   }

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
   * Assignment of LHS DeviceArray<Data> from RHS HostDArray<Data>.
   */
   template <typename Data>
   DeviceArray<Data>& 
   DeviceArray<Data>::operator = (const HostDArray<Data>& other)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("RHS HostDArray<Data> must be allocated.");
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
      // Make sure that owner array has not been deallocated / reallocated
      if (ownerPtr_) {
         UTIL_CHECK(ownerPtr_->isAllocated());
         UTIL_CHECK(ownerPtr_->capacity() == ownerCapacity_);
         UTIL_CHECK(ownerPtr_->cArray() == ownerDataPtr_);
      }

      return dataPtr_; 
   }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   const Data* DeviceArray<Data>::cArray() const
   {  
      // Make sure that owner array has not been deallocated / reallocated
      if (ownerPtr_) {
         UTIL_CHECK(ownerPtr_->isAllocated());
         UTIL_CHECK(ownerPtr_->capacity() == ownerCapacity_);
         UTIL_CHECK(ownerPtr_->cArray() == ownerDataPtr_);
      }

      return dataPtr_; 
   }

   /*
   * Return true if the array has been allocated, false otherwise.
   */
   template <typename Data>
   bool DeviceArray<Data>::isAllocated() const
   {  
      // Make sure that owner array has not been deallocated / reallocated
      if (ownerPtr_) {
         UTIL_CHECK(ownerPtr_->isAllocated());
         UTIL_CHECK(ownerPtr_->capacity() == ownerCapacity_);
         UTIL_CHECK(ownerPtr_->cArray() == ownerDataPtr_);
      }

      return (bool)dataPtr_; 
   }

   /*
   * Return allocated capacity.
   */
   template <typename Data>
   inline int DeviceArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Does this object own the underlying data? 
   */
   template <typename Data>
   inline bool DeviceArray<Data>::isOwner() const
   {  return ((nullptr == ownerPtr_) && (nullptr != dataPtr_)); }

} // namespace Pscf
#endif
