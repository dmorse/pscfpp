#ifndef PRDC_CUDA_R_FIELD_TPP
#define PRDC_CUDA_R_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"

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
    : DeviceArray<cudaReal>()
   {}

   /**
   * Allocating constructor.
   */
   template <int D>
   RField<D>::RField(IntVec<D> const & meshDimensions)
    : DeviceArray<cudaReal>()
   {  allocate(meshDimensions); }

   /*
   * Copy constructor.
   */
   template <int D>
   RField<D>::RField(const RField<D>& other)
    : DeviceArray<cudaReal>(other),
      meshDimensions_(0)
   {  meshDimensions_ = other.meshDimensions_; }

   /*
   * Destructor.
   */
   template <int D>
   RField<D>::~RField()
   {}

   /*
   * Allocate the underlying C array for data on a regular mesh.
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
      DeviceArray<cudaReal>::allocate(size);
   }

   /*
   * Associate this object with a slice of another DeviceArray.
   */
   template <int D>
   void RField<D>::associate(DeviceArray<cudaReal>& arr, int beginId, 
                             IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      DeviceArray<cudaReal>::associate(arr, beginId, size);
   }

   /*
   * Assignment from another RField<D>.
   */
   template <int D>
   RField<D>& RField<D>::operator = (const RField<D>& other)
   {
      DeviceArray<cudaReal>::operator = (other);
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
      DeviceArray<cudaReal>::operator = (other);

      return *this;
   }

}
}
}
#endif
