#ifndef PSCF_DEVICE_D_ARRAY_TPP
#define PSCF_DEVICE_D_ARRAY_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceDArray.h"
#include "HostDArray.h"
#include <cuda_runtime.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   DeviceDArray<Data>::DeviceDArray()
    : data_(0),
      capacity_(0)
   {}

   /*
   * Default constructor.
   */
   template <typename Data>
   DeviceDArray<Data>::DeviceDArray(int capacity)
    : data_(0),
      capacity_(0)
   {  allocate(capacity); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <typename Data>
   DeviceDArray<Data>::DeviceDArray(const DeviceDArray<Data>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other array must be allocated.");
      }

      allocate(other.capacity_);
      cudaMemcpy(data_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

   }

   /*
   * Destructor.
   */
   template <typename Data>
   DeviceDArray<Data>::~DeviceDArray()
   {
      if (isAllocated()) {
         cudaFree(data_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the array is already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <typename Data>
   void DeviceDArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate an array");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      gpuErrchk(cudaMalloc((void**) &data_, capacity * sizeof(Data)));
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this array is not allocated.
   */
   template <typename Data>
   void DeviceDArray<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      cudaFree(data_);
      data_ = 0; // reset to null ptr
      capacity_ = 0;
   }

   /*
   * Assignment from another DeviceDArray<Data>.
   */
   template <typename Data>
   DeviceDArray<Data>& 
   DeviceDArray<Data>::operator = (const DeviceDArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - RHS DeviceDArray<Data> must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("Other DeviceDArray<Data> must be allocated.");
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
      cudaMemcpy(data_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

      return *this;
   }

   /*
   * Assignment of LHS DeviceDArray<Data> from RHS HostDArray<Data>.
   */
   template <typename Data>
   DeviceDArray<Data>& DeviceDArray<Data>::operator = (const HostDArray<Data>& other)
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
      cudaMemcpy(data_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToDevice);

      return *this;
   }

} // namespace Pscf
#endif
