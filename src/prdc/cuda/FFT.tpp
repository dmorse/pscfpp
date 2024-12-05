#ifndef PRDC_CUDA_FFT_TPP
#define PRDC_CUDA_FFT_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"
#include <pscf/cuda/GpuResources.h>
#include <pscf/cuda/LinearAlgebra.h>

//forward declaration
//static __global__ void scaleRealData(cudaReal* data, rtype scale, int size);

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFT<D>::FFT()
    : meshDimensions_(0),
      rSize_(0),
      kSize_(0),
      rcfPlan_(0),
      criPlan_(0),
      ccPlan_(0),
      isSetup_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FFT<D>::~FFT()
   {
      if (rcfPlan_) {
         cufftDestroy(rcfPlan_);
      }
      if (criPlan_) {
         cufftDestroy(criPlan_);
      }
      if (ccPlan_) {
         cufftDestroy(ccPlan_);
      }
   }

   /*
   * Setup mesh dimensions.
   */
   template <int D>
   void FFT<D>::setup(IntVec<D> const & meshDimensions)
   {
      // Precondition
      UTIL_CHECK(!isSetup_);

      // Set mesh dimensions and sizes
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         rSize_ *= meshDimensions[i];
         if (i < D - 1) {
            kSize_ *= meshDimensions[i];
         } else {
            kSize_ *= (meshDimensions[i]/2 + 1);
         }
      }

      // Allocate rFieldCopy_ array if necessary
      if (!rFieldCopy_.isAllocated()) {
         rFieldCopy_.allocate(meshDimensions);
      } else {
         if (rFieldCopy_.capacity() != rSize_) {
            rFieldCopy_.deallocate();
            rFieldCopy_.allocate(meshDimensions);
         }
      }
      UTIL_CHECK(rFieldCopy_.capacity() == rSize_);

      // Allocate cFieldCopy_ array if necessary
      if (!cFieldCopy_.isAllocated()) {
         cFieldCopy_.allocate(meshDimensions);
      } else {
         if (cFieldCopy_.capacity() != rSize_) {
            cFieldCopy_.deallocate();
            cFieldCopy_.allocate(meshDimensions);
         }
      }
      UTIL_CHECK(cFieldCopy_.capacity() == rSize_);

      // Allocate kFieldCopy_ array if necessary
      if (!kFieldCopy_.isAllocated()) {
         kFieldCopy_.allocate(meshDimensions);
      } else {
         if (kFieldCopy_.capacity() != rSize_) {
            kFieldCopy_.deallocate();
            kFieldCopy_.allocate(meshDimensions);
         }
      }
      UTIL_CHECK(kFieldCopy_.capacity() == kSize_);

      // Make FFTW plans (explicit specializations)
      makePlans();

      isSetup_ = true;
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(RField<D> & rField, RFieldDft<D>& kField)
   const
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(rSize_, nBlocks, nThreads);

      // Check dimensions or setup
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);

      // Rescale outputted data. 
      cudaReal scale = 1.0/cudaReal(rSize_);
      //scaleRealData<<<nBlocks, nThreads>>>(rField.cArray(), scale, rSize_);
      scaleReal<<<nBlocks, nThreads>>>(rField.cArray(), scale, rSize_);
      
      // Perform transform
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecR2C(rcfPlan_, rField.cArray(), kField.cArray());
      #else
      result = cufftExecD2Z(rcfPlan_, rField.cArray(), kField.cArray());
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft real-to-complex forward transform");
      }
   }

   /*
   * Execute forward transform without destroying input.
   */
   template <int D>
   void FFT<D>::forwardTransformSafe(RField<D> const & rField, 
                                     RFieldDft<D>& kField)
   const
   {
      UTIL_CHECK(rFieldCopy_.capacity() == rField.capacity());

      rFieldCopy_ = rField;
      forwardTransform(rFieldCopy_, kField);
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFT<D>::inverseTransform(RFieldDft<D> & kField, RField<D>& rField) 
   const
   {
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2R(criPlan_, kField.cArray(), rField.cArray());
      #else
      result = cufftExecZ2D(criPlan_, kField.cArray(), rField.cArray());
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW( "Failure in cufft complex-to-real inverse transform");
      }
   
   }

   /*
   * Execute inverse (complex-to-real) transform without destroying input.
   */
   template <int D>
   void FFT<D>::inverseTransformSafe(RFieldDft<D> const & kField, 
                                     RField<D>& rField) const
   {
      UTIL_CHECK(kFieldCopy_.capacity()==kField.capacity());

      kFieldCopy_ = kField;
      inverseTransform(kFieldCopy_, rField);
   }

   // Complex-to-Complex Transforms

   /*
   * Execute forward complex-to-complex transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(CField<D> & rField, CField<D>& kField)
   const
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(rSize_, nBlocks, nThreads);

      // Check dimensions or setup
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);

      // Rescale outputted data. 
      cudaReal scale = 1.0/cudaReal(rSize_);
      //scaleComplexData<<<nBlocks, nThreads>>>(rField.cArray(), scale, rSize_);
      scaleComplex<<<nBlocks, nThreads>>>(rField.cArray(), scale, rSize_);
      
      // Perform transform
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2C(ccPlan_, rField.cArray(), kField.cArray(), CUFFT_FORWARD);
      #else
      result = cufftExecZ2Z(ccPlan_, rField.cArray(), kField.cArray(), CUFFT_FORWARD);
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft real-to-complex forward transform");
      }
   }

   /*
   * Execute forward transform without destroying input.
   */
   template <int D>
   void FFT<D>::forwardTransformSafe(CField<D> const & rField, 
                                     CField<D>& kField)
   const
   {
      UTIL_CHECK(cFieldCopy_.capacity() == rField.capacity());

      cFieldCopy_ = rField;
      forwardTransform(cFieldCopy_, kField);
   }

   /*
   * Execute inverse (complex-to-complex) transform.
   */
   template <int D>
   void FFT<D>::inverseTransform(CField<D> & kField, CField<D>& rField) 
   const
   {
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2C(ccPlan_, kField.cArray(), rField.cArray(), CUFFT_INVERSE);
      #else
      result = cufftExecZ2Z(ccPlan_, kField.cArray(), rField.cArray(), CUFFT_INVERSE);
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW( "Failure in cufft complex-to-real inverse transform");
      }
   
   }

   /*
   * Execute inverse (complex-to-complex) transform without destroying input.
   */
   template <int D>
   void FFT<D>::inverseTransformSafe(CField<D> const & kField, 
                                     CField<D>& rField) const
   {
      UTIL_CHECK(cFieldCopy_.capacity() == kField.capacity());

      cFieldCopy_ = kField;
      inverseTransform(cFieldCopy_, rField);
   }

}
}
}
#endif
