#ifndef PRDC_CUDA_C_FIELD_TPP
#define PRDC_CUDA_C_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CField.h"
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
   CField<D>::CField()
    : DeviceDArray<cudaComplex>()
   {}

   /**
   * Allocating constructor.
   */
   template <int D>
   CField<D>::CField(IntVec<D> const & meshDimensions)
    : DeviceDArray<cudaComplex>()
   {  allocate(meshDimensions); }

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
   *
   *\param other the Field to be copied.
   */
   template <int D>
   CField<D>::CField(const CField<D>& other)
    : DeviceDArray<cudaComplex>(other),
      meshDimensions_(0)
   {
      meshDimensions_ = other.meshDimensions_;
   }

   /*
   * Assignment from another RField<D>.
   */
   template <int D>
   CField<D>& CField<D>::operator = (const CField<D>& other)
   {
      DeviceDArray<cudaComplex>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

   /*
   * Assignment from RHS HostDArray<Data> host array.
   */
   template <int D>
   CField<D>& CField<D>::operator = (HostDArray<cudaComplex> const & other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostDArray<cudaComplex> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS CField<D> is not allocated.");
      }
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Use base class assignment operator to copy elements
      DeviceDArray<cudaComplex>::operator = (other);

      return *this;
   }

   /*
   * Allocate the underlying C array sized for an associated mesh.
   */
   template <int D>
   void CField<D>::allocate(IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      DeviceDArray<cudaComplex>::allocate(size);
   }

}
}
}
#endif
