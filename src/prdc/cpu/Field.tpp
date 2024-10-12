#ifndef PRDC_CPU_FIELD_TPP
#define PRDC_CPU_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
#include <util/misc/Memory.h>

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

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
   * Destructor.
   */
   template <typename Data>
   Field<Data>::~Field()
   {
      if (isAllocated()) {
         fftw_free(data_);
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
         UTIL_THROW("Attempt to allocate Field with capacity <= 0");
      }
      capacity_ = capacity;
      data_ = (Data*) fftw_malloc(sizeof(Data)*capacity);
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
      fftw_free(data_);
      capacity_ = 0;
   }

   #ifdef PSCF_CUDA
   /*
   * Assignment of Cpu::Field<Data> LHS from Cuda:RField<Data> RHS.
   *
   * This is a member of Cpu::Field<Data> that is declared in the header for that
   * class template, but defined here, and so is only compiled or used when 
   * compilation of CUDA code is enabled. 
   *
   */
   template <typename Data>
   Field<Data>& Field<Data>::operator = (const Cuda::Field<Data>& other)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("RHS Cuda::Field<D> must be allocated.");
      }

      // Allocate this if necessary, otherwise require equal dimensions
      if (!isAllocated()) {
         allocate(other.capacity());
      } else if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Copy elements
      cudaMemcpy(data_, other.cField(), 
                 capacity_ * sizeof(Data), cudaMemcpyDeviceToHost);

      return *this;
   }
   #endif

}
}
}
#endif
