#ifndef PSPG_DFIELD_TPP
#define PSPG_DFIELD_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DField.h"
#include <pspg/math/GpuResources.h>
#include <cuda_runtime.h>



namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   DField<Data>::DField()
    : data_(0),
      capacity_(0)
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   DField<Data>::~DField()
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
   void DField<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a DField");
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
   void DField<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      cudaFree(data_);
      capacity_ = 0;
   }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <typename Data>
   DField<Data>::DField(const DField<Data>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }

      allocate(other.capacity_);
      cudaMemcpy(data_, other.cDField(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

   }

   /*
   * Assignment.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other Field is not allocated.
   * \throw Exception if both Fields are allocated with unequal capacities.
   *
   * \param other the rhs Field
   */
   template <typename Data>
   DField<Data>& DField<Data>::operator = (const DField<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }

      if (!isAllocated()) {
         allocate(other.capacity());
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(data_, other.cDField(), capacity_ * sizeof(Data), cudaMemcpyDeviceToDevice);

      return *this;
   }

}
}
#endif
