#ifndef PRDC_CUDA_R_FIELD_TPP
#define PRDC_CUDA_R_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include <pscf/cuda/HostDArray.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf;

   /**
   * Default constructor.
   */
   template <int D>
   RField<D>::RField()
    : DeviceDArray<cudaReal>()
   {}

   /**
   * Allocating constructor.
   */
   template <int D>
   RField<D>::RField(IntVec<D> const & meshDimensions)
    : DeviceDArray<cudaReal>()
   {  allocate(meshDimensions); }

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <int D>
   RField<D>::RField(const RField<D>& other)
    : DeviceDArray<cudaReal>(other),
      meshDimensions_(0)
   {  meshDimensions_ = other.meshDimensions_; }

   /*
   * Destructor.
   */
   template <int D>
   RField<D>::~RField()
   {}

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RField<D>::allocate(IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      DeviceDArray<cudaReal>::allocate(size);
   }

   /*
   * Assignment from another RField<D>.
   */
   template <int D>
   RField<D>& RField<D>::operator = (const RField<D>& other)
   {
      DeviceDArray<cudaReal>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

   /*
   * Assignment of this RField<D> from RHS HostDArray<Data> host array.
   */
   template <int D>
   RField<D>& RField<D>::operator = (const HostDArray<cudaReal>& other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostDArray<cudaReal> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS RField<D> is not allocated.");
      }
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Use base class assignment operator to copy elements
      DeviceDArray<cudaReal>::operator = (other);

      return *this;
   }

}
}
}
#endif
