#ifndef PSCF_HOST_D_ARRAY_TPP
#define PSCF_HOST_D_ARRAY_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostDArray.h"
#include "DeviceArray.h"

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   HostDArray<Data>::HostDArray()
    : dataPtr_(nullptr),
      capacity_(0)
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   HostDArray<Data>::HostDArray(int capacity)
    : dataPtr_(nullptr),
      capacity_(0)
   {  allocate(capacity); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the HostDArray to be copied.
   */
   template <typename Data>
   HostDArray<Data>::HostDArray(const HostDArray<Data>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other HostDArray must be allocated.");
      }

      allocate(other.capacity_);
      cudaMemcpy(dataPtr_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToHost);

   }

   /*
   * Destructor.
   */
   template <typename Data>
   HostDArray<Data>::~HostDArray()
   {
      if (isAllocated()) {
         cudaFreeHost(dataPtr_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   template <typename Data>
   void HostDArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a HostDArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      gpuErrChk(cudaMallocHost((void**) &dataPtr_, capacity * sizeof(Data)));
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this HostDArray is not allocated.
   */
   template <typename Data>
   void HostDArray<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Attempt to deallocate unallocated HostDArray");
      }
      cudaFreeHost(dataPtr_);
      dataPtr_ = nullptr; // reset to null
      capacity_ = 0;
   }

   /*
   * Assignment from another HostDArray<Data> host array.
   */
   template <typename Data>
   HostDArray<Data>& 
   HostDArray<Data>::operator = (const HostDArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - RHS array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("Other HostDArray must be allocated.");
      }

      // Allocate this LHS array if necessary 
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacities
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign HostDArrays of unequal capacity");
      }

      // Copy elements from RHS to LHS
      cudaMemcpy(dataPtr_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToHost);

      return *this;
   }

   /*
   * Assignment from a DeviceArray<Data> RHS device array.
   */
   template <typename Data>
   HostDArray<Data>& 
   HostDArray<Data>::operator = (const DeviceArray<Data>& other)
   {
      // Precondition - RHS array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DeviceArray<Data> must be allocated.");
      }

      // Allocate this if necessary 
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacities
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy all elements
      cudaMemcpy(dataPtr_, other.cArray(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToHost);

      return *this;
   }

}
#endif
