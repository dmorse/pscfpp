#ifndef PSCF_DEVICE_ARRAY_H
#define PSCF_DEVICE_ARRAY_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include <util/global.h>
#include <cufft.h>

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
      * Does this object own the underlying array? 
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
      Data* data_;

      /// Allocated size (capacity) of the data_ array.
      int capacity_;

      /// Does this object own the underlying array? 
      bool isOwner_;

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
      Data const * ownerData_;
   };

   /*
   * Return allocated capacity.
   */
   template <typename Data>
   inline int DeviceArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Does this object own the underlying array? 
   */
   template <typename Data>
   inline bool DeviceArray<Data>::isOwner() const
   {  return isOwner_; }

   #ifndef PSCF_DEVICE_ARRAY_TPP
   extern template class DeviceArray<cudaReal>;
   extern template class DeviceArray<cudaComplex>;
   #endif

} // namespace Pscf
#endif
