#ifndef PRDC_CPU_R_FIELD_TPP
#define PRDC_CPU_R_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include "Field.tpp"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RField<D>::RField()
    : Field<double>(),
      meshDimensions_(0)
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
      if (other.isAllocated()) {
         allocate(other.meshDimensions_);
         for (int i = 0; i < capacity_; ++i) {
            data_[i] = other.data_[i];
         }
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
         allocate(other.meshDimensions_);
      } 
      UTIL_CHECK(capacity_ == other.capacity_);
      UTIL_CHECK(meshDimensions_ == other.meshDimensions_);

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RField<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      Field<double>::allocate(size);
   }

   /*
   * Dellocate the underlying C array and clear meshDimensions.
   */
   template <int D>
   void RField<D>::deallocate()
   {
      Field<double>::deallocate();
      for (int i = 0; i < D; ++i) {
         meshDimensions_[i] = 0;
      }
   }

}
}
}
#endif
