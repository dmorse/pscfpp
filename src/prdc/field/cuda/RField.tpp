#ifndef PRDC_CUDA_R_FIELD_TPP
#define PRDC_CUDA_R_FIELD_TPP

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf;

   /**
   * Default constructor.
   */
   template <int D>
   RField<D,CUT>::RField()
    : DeviceArray<cudaReal>()
   {}

   /**
   * Allocating constructor.
   */
   template <int D>
   RField<D,CUT>::RField(IntVec<D> const & meshDimensions)
    : DeviceArray<cudaReal>()
   {  allocate(meshDimensions); }

   /*
   * Copy constructor.
   */
   template <int D>
   RField<D,CUT>::RField(RField<D,CUT> const & other)
    : DeviceArray<cudaReal>(other),
      meshDimensions_(0)
   {  meshDimensions_ = other.meshDimensions_; }

   /*
   * Destructor.
   */
   template <int D>
   RField<D,CUT>::~RField()
   {}

   /*
   * Allocate the underlying C array for data on a regular mesh.
   */
   template <int D>
   void RField<D,CUT>::allocate(IntVec<D> const & meshDimensions)
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
   void RField<D,CUT>::associate(DeviceArray<cudaReal>& arr, int beginId, 
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
   * Assignment from another RField.
   */
   template <int D>
   RField<D,CUT>& 
   RField<D,CUT>::operator = (const RField<D,CUT>& other)
   {
      DeviceArray<cudaReal>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

   /*
   * Assignment of this RField from RHS HostDArray.
   */
   template <int D>
   RField<D,CUT>& 
   RField<D,CUT>::operator = (const HostDArray<cudaReal>& other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostDArray<cudaReal> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS RField is not allocated.");
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
#endif
