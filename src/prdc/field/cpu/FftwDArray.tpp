#ifndef PRDC_CPU_FFTW_D_ARRAY_TPP
#define PRDC_CPU_FFTW_D_ARRAY_TPP

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FftwDArray.h"
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
         try {
            fftw_free(data_);
            Memory::sub<Data>(capacity_);
         } catch (...) {
            std::cout << "Exception in FftwDArray destructor";
         }
      }
      data_ = nullptr;
      capacity_ = 0;
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
      data_ = (Data*) fftw_malloc(sizeof(Data)*capacity);
      capacity_ = capacity;
      Memory::add<Data>(capacity);
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
      data_ = nullptr;
      Memory::sub<Data>(capacity_);
      capacity_ = 0;
   }

}
}
}
#endif
