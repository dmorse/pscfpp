#ifndef PRDC_CPU_FFTW_D_ARRAY_TPP
#define PRDC_CPU_FFTW_D_ARRAY_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FftwDArray.h"
#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   FftwDArray<Data>::FftwDArray()
    : Array<Data>()
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   FftwDArray<Data>::~FftwDArray()
   {
      if (isAllocated()) {
         fftw_free(data_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the FftwDArray has already been allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <typename Data>
   void FftwDArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a FftwDArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate FftwDArray with capacity <= 0");
      }
      capacity_ = capacity;
      data_ = (Data*) fftw_malloc(sizeof(Data)*capacity);
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this FftwDArray is not allocated.
   */
   template <typename Data>
   void FftwDArray<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      fftw_free(data_);
      capacity_ = 0;
   }

}
}
}
#endif
