/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include <util/misc/Memory.h>

#include <fftw3.h>

namespace Pssp
{

   using namespace Util;

   /*
   * Default constructor.
   */
   RField::RField()
    : data_(0),
      capacity_(0)
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the RField to be copied.
   */
   RField::RField(const RField& other)
    : data_(0),
      capacity_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other RField must be allocated.");
      }
      data_ = (double*) fftw_malloc(sizeof(double)*other.capacity_);
      //Memory::allocate(data_, other.capacity_);
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Destructor.
   */
   RField::~RField()
   {
      if (isAllocated()) {
         fftw_free(data_);
         // Memory::deallocate<double>(data_, capacity_);
         capacity_ = 0;
      }
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other RField is not allocated.
   * \throw Exception if both RFields are allocated with unequal capacities.
   *
   * \param other the rhs RField
   */
   RField& RField::operator = (const RField& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other RField must be allocated.");
      }

      if (!isAllocated()) {
         allocate(other.capacity());
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign RFields of unequal capacity");
      }

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the RField has already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   void RField::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a RField");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      data_ = (double*) fftw_malloc(sizeof(double)*capacity);
      //Memory::allocate<double>(data_, capacity);
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this RField is not allocated.
   */
   void RField::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      fftw_free(data_);
      //Memory::deallocate<double>(data_, capacity_);
      capacity_ = 0;
   }

}
