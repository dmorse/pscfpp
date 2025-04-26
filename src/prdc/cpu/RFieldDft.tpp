#ifndef PRDC_CPU_R_FIELD_DFT_TPP
#define PRDC_CPU_R_FIELD_DFT_TPP

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.h"
#include <prdc/cpu/FFT.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   RFieldDft<D>::RFieldDft()
    : FftwDArray<fftw_complex>(),
      meshDimensions_(0),
      dftDimensions_(0)
   {}

   /*
   * Destructor.
   */
   template <int D>
   RFieldDft<D>::~RFieldDft()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <int D>
   RFieldDft<D>::RFieldDft(const RFieldDft<D>& other)
    : FftwDArray<fftw_complex>(),
      meshDimensions_(0),
      dftDimensions_(0)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field must be allocated.");
      }
      allocate(other.meshDimensions_);
      for (int i = 0; i < capacity_; ++i) {
         data_[i][0] = other.data_[i][0];
         data_[i][1] = other.data_[i][1];
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
   RFieldDft<D>& RFieldDft<D>::operator = (const RFieldDft<D>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other RFieldDft must be allocated.");
      }

      // Allocate if necessary, check dimensions
      if (!isAllocated()) {
         allocate(other.meshDimensions_);
      }
      UTIL_CHECK(capacity_ == other.capacity_);
      UTIL_CHECK(meshDimensions_ == other.meshDimensions_);
      UTIL_CHECK(dftDimensions_ == other.dftDimensions_);

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i][0] = other.data_[i][0];
         data_[i][1] = other.data_[i][1];
      }

      return *this;
   }

   /*
   * Allocate the underlying array for an FFT grid.
   */
   template <int D>
   void RFieldDft<D>::allocate(IntVec<D> const & meshDimensions)
   {
      // Copy real space grid dimensions
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
      }

      // Compute dimensions and size for k-space grid
      int size; 
      FFT<D>::computeKMesh(meshDimensions, dftDimensions_, size);

      // Allocate memory
      FftwDArray<fftw_complex>::allocate(size);
   }

   /*
   * Dellocate the underlying C array and clear dimensions.
   */
   template <int D>
   void RFieldDft<D>::deallocate()
   {
      FftwDArray<fftw_complex>::deallocate();
      for (int i = 0; i < D; ++i) {
         meshDimensions_[i] = 0;
         dftDimensions_[i] = 0;
      }
   }

}
}
}
#endif
