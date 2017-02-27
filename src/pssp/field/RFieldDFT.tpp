#ifndef PSSP_R_FIELD_DFT_TPP
#define PSSP_R_FIELD_DFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDFT.h"

namespace Pssp
{

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RFieldDFT<D>::RFieldDFT()
    : Field<fftw_complex>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   RFieldDFT<D>::~RFieldDFT()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the RField<D> to be copied.
   */
   template <int D>
   RFieldDFT<D>::RFieldDFT(const RFieldDFT<D>& other)
    : Field<fftw_complex>()
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }
      data_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*other.capacity_);
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i][0] = other.data_[i][0];
         data_[i][1] = other.data_[i][1];
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
   RFieldDFT<D>& RFieldDFT<D>::operator = (const RFieldDFT<D>& other)
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
         data_[i][0] = other.data_[i][0];
         data_[i][1] = other.data_[i][1];
      }
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

}
#endif
