#ifndef PRDC_CPU_FIELD_TPP
#define PRDC_CPU_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
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
    : Array<Data>()
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
   * Throw an Exception if the Field has already been allocated.
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

}
}
}
#endif
