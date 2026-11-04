#ifndef PSCF_R_FIELD_CU_TPP
#define PSCF_R_FIELD_CU_TPP

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include <pscf/backend/cuda/HostSendArray.h>

namespace Pscf {

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RField<D,CUT>::RField()
    : BackendArray<cudaReal,CUT>()
   {}

   /**
   * Allocating constructor.
   */
   template <int D>
   RField<D,CUT>::RField(IntVec<D> const & meshDimensions)
    : BackendArray<cudaReal,CUT>()
   {  allocate(meshDimensions); }

   /*
   * Copy constructor.
   */
   template <int D>
   RField<D,CUT>::RField(RField<D,CUT> const & other)
    : BackendArray<cudaReal,CUT>(other),
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
      BackendArray<cudaReal,CUT>::allocate(size);
   }

   /*
   * Deallocate the underlying C array and zero grid dimensions.
   */
   template <int D>
   void RField<D,CUT>::deallocate()
   {
      BackendArray<cudaReal,CUT>::deallocate();
      for (int i = 0; i < D; ++i) {
         meshDimensions_[i] = 0;
      }
   }

   /*
   * Associate this object with a slice of another BackendArray.
   */
   template <int D>
   void RField<D,CUT>::associate(BackendArray<cudaReal,CUT>& arr, 
                                 int beginId, 
                                 IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      BackendArray<cudaReal,CUT>::associate(arr, beginId, size);
   }

   /*
   * Assignment from another RField.
   */
   template <int D>
   RField<D,CUT>& 
   RField<D,CUT>::operator = (const RField<D,CUT>& other)
   {
      BackendArray<cudaReal,CUT>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

   /*
   * Assignment of this RField from RHS HostSendArray.
   */
   template <int D>
   RField<D,CUT>& 
   RField<D,CUT>::operator = (const HostSendArray<cudaReal,CUT>& other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      UTIL_CHECK(other.isAllocated());
      UTIL_CHECK((BackendArray<cudaReal,CUT>::isAllocated()));
      UTIL_CHECK((BackendArray<cudaReal,CUT>::capacity() == other.capacity()));

      // Use base class assignment operator to copy elements
      BackendArray<cudaReal,CUT>::operator = (other);

      return *this;
   }

}
#endif
