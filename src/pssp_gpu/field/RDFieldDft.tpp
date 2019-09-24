#ifndef PSSP_GPU_R_DFIELD_DFT_TPP
#define PSSP_GPU_R_DFIELD_DFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RDFieldDft.h"

namespace Pscf {
namespace Pssp_gpu {

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RDFieldDft<D>::RDFieldDft()
    : DField<cufftComplex>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   RDFieldDft<D>::~RDFieldDft()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the RField<D> to be copied.
   */
   template <int D>
   RDFieldDft<D>::RDFieldDft(const RDFieldDft<D>& other)
    : DField<cufftComplex>()
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }

      capacity_ = other.capacity_;
      cudaMalloc((void**) &data_, capacity_ * sizeof(cufftComplex));
      
      cudaMemcpy(data_, other.cDField(), capacity_ * sizeof(cufftComplex), cudaMemcpyDeviceToDevice);
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
   RDFieldDft<D>& RDFieldDft<D>::operator = (const RDFieldDft<D>& other)
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
      cudaMemcpy(data_, other.cDField(), capacity_ * sizeof(cufftComplex), cudaMemcpyDeviceToDevice);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

}
}
#endif
