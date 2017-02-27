#ifndef PSSP_R_FIELD_TPP
#define PSSP_R_FIELD_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"

namespace Pssp
{

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RField<D>::RField()
    : Field<double>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   RField<D>::~RField()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <int D>
   RField<D>::RField(const RField<D>& other)
    : Field<double>(),
      meshDimensions_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }
      data_ = (double*) fftw_malloc(sizeof(double)*other.capacity_);
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
      meshDimensions_ = other.meshDimensions_;
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
   template <int D>
   RField<D>& RField<D>::operator = (const RField<D>& other)
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
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

}
#endif
