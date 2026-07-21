#ifndef PRDC_R_FIELD_DFT_CU_TPP
#define PRDC_R_FIELD_DFT_CU_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.h"
#include "FFT.h"

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   RFieldDft<D, CudaTp<D> >::RFieldDft()
    : DeviceArray<cudaComplex>()
   {}

   /*
   * Allocating constructor (calls allocate).
   */
   template <int D>
   RFieldDft<D, CudaTp<D> >::RFieldDft(IntVec<D> const & meshDimensions)
    : DeviceArray<cudaComplex>()
   {  allocate(meshDimensions); }

   /*
   * Destructor.
   */
   template <int D>
   RFieldDft<D, CudaTp<D> >::~RFieldDft()
   {}

   /*
   * Copy constructor.
   */
   template <int D>
   RFieldDft<D, CudaTp<D> >::RFieldDft(const RFieldDft<D, CudaTp<D> >& other)
    : DeviceArray<cudaComplex>(other)
   {
      meshDimensions_ = other.meshDimensions_;
      dftDimensions_ = other.dftDimensions_;
   }

   /*
   * Assignment from another RFieldDft.
   */
   template <int D>
   RFieldDft<D, CudaTp<D> >&
   RFieldDft<D, CudaTp<D> >::operator = (RFieldDft<D, CudaTp<D> > const & other)
   {
      // Assign data and size of underlying array
      DeviceArray<cudaComplex>::operator = (other);

      // Assign more specialized data members
      meshDimensions_ = other.meshDimensions_;
      dftDimensions_ = other.dftDimensions_;

      return *this;
   }

   /*
   * Assignment from RHS HostDArray<cudaComplex>.
   */
   template <int D>
   RFieldDft<D, CudaTp<D> >&
   RFieldDft<D, CudaTp<D> >::operator = (HostDArray<cudaComplex> const & other)
   {
      // Preconditions: Both arrays must be allocated with equal capacities
      if (!other.isAllocated()) {
         UTIL_THROW("Error: RHS HostDArray<cudaComplex> is not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("Error: LHS RFieldDft<D, CudaTp<D> > is not allocated.");
      }
      if (capacity_ != other.capacity()) {
         UTIL_THROW("Cannot assign Fields of unequal capacity");
      }

      // Use base class assignment operator to copy elements
      DeviceArray<cudaComplex>::operator = (other);

      return *this;
   }

   /*
   * Allocate underlying DeviceArray<cudaComplex> for the DFT mesh.
   */
   template <int D>
   void RFieldDft<D, CudaTp<D> >::allocate(const IntVec<D>& meshDimensions)
   {
      // Copy and validate dimensions of real space grid
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
      }

      // Compute dimensions and size of Fourier space mesh
      int size;
      FFT<D, CudaTp<D> >::computeKMesh(meshDimensions, dftDimensions_, size);

      // Allocate complex array on the GPU with size of DFT mesh
      DeviceArray<cudaComplex>::allocate(size);
   }

   /*
   * Associate this with a slice of a DeviceArray<cudaComplex>.
   */
   template <int D>
   void RFieldDft<D, CudaTp<D> >::associate(
                DeviceArray<cudaComplex>& arr,
                int beginId,
                IntVec<D> const & meshDimensions)
   {
      // Copy and validate dimensions of real space grid
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
      }

      // Compute dimensions and size of Fourier space mesh
      int size;
      FFT<D, CudaTp<D> >::computeKMesh(meshDimensions, dftDimensions_, size);

      // Associate data with a slice of input array arr
      DeviceArray<cudaComplex>::associate(arr, beginId, size);
   }

} // namespace Prdc
} // namespace Pscf
#endif
