#ifndef PRDC_CUDA_R_FIELD_DFT_TPP
#define PRDC_CUDA_R_FIELD_DFT_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.h"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   RFieldDft<D>::RFieldDft()
    : DeviceArray<cudaComplex>()
   {}

   /*
   * Allocating constructor (calls allocate).
   */
   template <int D>
   RFieldDft<D>::RFieldDft(IntVec<D> const & meshDimensions)
    : DeviceArray<cudaComplex>()
   {  allocate(meshDimensions); }

   /*
   * Destructor.
   */
   template <int D>
   RFieldDft<D>::~RFieldDft()
   {}

   /*
   * Copy constructor.
   */
   template <int D>
   RFieldDft<D>::RFieldDft(const RFieldDft<D>& other)
    : DeviceArray<cudaComplex>(other)
   {
      meshDimensions_ = other.meshDimensions_;
      dftDimensions_ = other.dftDimensions_;
   }

   /*
   * Assignment, element-by-element from another RFieldDft<D>.
   */
   template <int D>
   RFieldDft<D>& RFieldDft<D>::operator = (const RFieldDft<D>& other)
   {
      
      DeviceArray<cudaComplex>::operator = (other);
      meshDimensions_ = other.meshDimensions_;
      dftDimensions_ = other.dftDimensions_;

      return *this;
   }

   /*
   * Assignment of RFieldDft<D> from RHS HostDArray<Data> host array.
   */
   template <int D>
   RFieldDft<D>& 
   RFieldDft<D>::operator = (const HostDArray<cudaComplex>& other)
   {
      // Preconditions: both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostDArray<cudaComplex> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS RFieldDft<D> is not allocated.");
      }
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Use base class assignment operator to copy elements
      DeviceArray<cudaComplex>::operator = (other);

      return *this;
   }

   /*
   * Allocate the underlying C array for the dftDimensions mesh.
   */
   template <int D>
   void RFieldDft<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         if (i < D - 1) {
            dftDimensions_[i] = meshDimensions[i];
            size *= meshDimensions[i];
         } else {
            dftDimensions_[i] = (meshDimensions[i]/2 + 1); 
            size *= (meshDimensions[i]/2 + 1);
         }
      }
      // Note: size denotes size of mesh with dftDimensions 
      DeviceArray<cudaComplex>::allocate(size);
   }

   /*
   * Associate this object with a slice of another DeviceArray.
   */
   template <int D>
   void RFieldDft<D>::associate(DeviceArray<cudaComplex>& arr, int beginId, 
                                IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         if (i < D - 1) {
            dftDimensions_[i] = meshDimensions[i];
            size *= meshDimensions[i];
         } else {
            dftDimensions_[i] = (meshDimensions[i]/2 + 1); 
            size *= (meshDimensions[i]/2 + 1);
         }
      }
      // Note: size denotes size of mesh with dftDimensions
      DeviceArray<cudaComplex>::associate(arr, beginId, size);
   }

}
}
}
#endif
