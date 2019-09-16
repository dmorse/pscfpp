#ifndef PSSP_GPU_DFIELD_TPP
#define PSSP_GPU_DFIELD_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DField.h"
#include <cuda_runtime.h>



namespace Pscf {
namespace Pssp_gpu
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
      cudaMalloc((void**) &data_, capacity * sizeof(Data));
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

}
}
#endif
