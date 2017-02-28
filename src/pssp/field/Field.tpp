#ifndef PSSP_FIELD_TPP
#define PSSP_FIELD_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
#include <util/misc/Memory.h>

#include <fftw3.h>

namespace Pscf {
namespace Pssp
{

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

   #if 0
   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <typename Data>
   Field<Data>::Field(const Field<Data>& other)
    : data_(0),
      capacity_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }
      data_ = (Data*) fftw_malloc(sizeof(Data)*other.capacity_);
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other Field is not allocated.
   * \throw Exception if both Fields are allocated with unequal capacities.
   *
   * \param other the rhs Field
   */
   template <typename Data>
   Field<Data>& Field<Data>::operator = (const Field<Data>& other)
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
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }
   #endif

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
      data_ = (Data*) fftw_malloc(sizeof(Data)*capacity);
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
      fftw_free(data_);
      capacity_ = 0;
   }

}
}
#endif
