#ifndef PRDC_CUDA_FIELD_TPP
#define PRDC_CUDA_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/Field.h>
#include <prdc/cuda/HostField.h>
#include <pscf/cuda/GpuResources.h>
#include <cuda_runtime.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   Field<Data>::Field()
    : data_(0),
      capacity_(0)
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <typename Data>
   Field<Data>::Field(const Field<Data>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }

      allocate(other.capacity_);
      cudaMemcpy(data_, other.cField(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

   }

   /*
   * Destructor.
   */
   template <typename Data>
   Field<Data>::~Field()
   {
      if (isAllocated()) {
         cudaFree(data_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the Field has already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <typename Data>
   void Field<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a Field");
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
   * Throw an Exception if this Field is not allocated.
   */
   template <typename Data>
   void Field<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      cudaFree(data_);
      capacity_ = 0;
   }

   /*
   * Assignment from another Cuda::Field<Data>.
   */
   template <typename Data>
   Field<Data>& Field<Data>::operator = (const Field<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - RHS Field<Data> must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field<Data> must be allocated.");
      }

      // If this is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacity values 
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(data_, other.cField(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

      return *this;
   }

   /*
   * Assignment of LHS Field<Data> from LHS HostField<Data> host array.
   */
   template <typename Data>
   Field<Data>& Field<Data>::operator = (const HostField<Data>& other)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("RHS HostField<Data> must be allocated.");
      }

      // Allocate this if necessary
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacity values
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(data_, other.cField(), 
                 capacity_ * sizeof(Data), cudaMemcpyHostToDevice);

      return *this;
   }

} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif
