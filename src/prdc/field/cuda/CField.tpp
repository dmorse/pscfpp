#ifndef PRDC_C_FIELD_CU_TPP
#define PRDC_C_FIELD_CU_TPP

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CField.h"

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf;

   /**
   * Default constructor.
   */
   template <int D>
   CField<D,CUT>::CField()
    : DeviceArray<cudaComplex,CUT>()
   {}

   /**
   * Allocating constructor.
   */
   template <int D>
   CField<D,CUT>::CField(IntVec<D> const & meshDimensions)
    : DeviceArray<cudaComplex,CUT>()
   {  allocate(meshDimensions); }

   /*
   * Destructor.
   */
   template <int D>
   CField<D,CUT>::~CField()
   {}

   /*
   * Copy constructor.
   */
   template <int D>
   CField<D,CUT>::CField(const CField<D,CUT>& other)
    : DeviceArray<cudaComplex,CUT>(other),
      meshDimensions_(0)
   {
      meshDimensions_ = other.meshDimensions_;
   }

   /*
   * Assignment from another RField<D,CUT>.
   */
   template <int D>
   CField<D,CUT>& CField<D,CUT>::operator = (const CField<D,CUT>& other)
   {
      DeviceArray<cudaComplex,CUT>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

   /*
   * Assignment from RHS HostArray<Data> host array.
   */
   template <int D>
   CField<D,CUT>& CField<D,CUT>::operator = (HostArray<cudaComplex,CUT> const & other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostArray<cudaComplex,CUT> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS CField<D,CUT> is not allocated.");
      }
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Use base class assignment operator to copy elements
      DeviceArray<cudaComplex,CUT>::operator = (other);

      return *this;
   }

   /*
   * Allocate the underlying C array sized for an associated mesh.
   */
   template <int D>
   void CField<D,CUT>::allocate(IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      DeviceArray<cudaComplex,CUT>::allocate(size);
   }

   /*
   * Associate this object with a slice of another DeviceArray.
   */
   template <int D>
   void CField<D,CUT>::associate(
                             DeviceArray<cudaComplex,CUT>& arr, 
                             int beginId, 
                             IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      DeviceArray<cudaComplex,CUT>::associate(arr, beginId, size);
   }

}
}
#endif
