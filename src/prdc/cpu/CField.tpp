#ifndef PRDC_CPU_C_FIELD_TPP
#define PRDC_CPU_C_FIELD_TPP

/*
* PSCF Package
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CField.h"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   CField<D>::CField()
    : FftwDArray<fftw_complex>(),
      meshDimensions_()
   {}

   /*
   * Destructor.
   */
   template <int D>
   CField<D>::~CField()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   */
   template <int D>
   CField<D>::CField(const CField<D>& other)
    : FftwDArray<fftw_complex>(),
      meshDimensions_()
   {
      if (other.isAllocated() && other.capacity_ > 0) {
         FftwDArray<fftw_complex>::allocate(other.capacity_);
         meshDimensions_ = other.meshDimensions_;
         for (int i = 0; i < capacity_; ++i) {
            data_[i][0] = other.data_[i][0];
            data_[i][1] = other.data_[i][1];
         }
      }
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   */
   template <int D>
   CField<D>& CField<D>::operator = (CField<D> const & other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other CField must be allocated in assignment.");
      }

      if (!isAllocated()) {
         allocate(other.meshDimensions_);
      }
      UTIL_CHECK(capacity_ == other.capacity_);
      UTIL_CHECK(meshDimensions_ == other.meshDimensions_);

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i][0] = other.data_[i][0];
         data_[i][1] = other.data_[i][1];
      }

      return *this;
   }

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void CField<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      FftwDArray<fftw_complex>::allocate(size);
   }

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void CField<D>::deallocate()
   {
      FftwDArray<fftw_complex>::deallocate();
      for (int i = 0; i < D; ++i) {
         meshDimensions_[i] = 0;
      }
   }

}
}
}
#endif
