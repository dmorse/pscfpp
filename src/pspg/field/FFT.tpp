#ifndef PSPG_FFT_TPP
#define PSPG_FFT_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"
#include <pspg/math/GpuResources.h>

//forward declaration
//static __global__ void scaleRealData(cudaReal* data, rtype scale, int size);

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFT<D>::FFT()
    : meshDimensions_(0),
      rSize_(0),
      kSize_(0),
      fPlan_(0),
      iPlan_(0),
      isSetup_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FFT<D>::~FFT()
   {
      if (fPlan_) {
         cufftDestroy(fPlan_);
      }
      if (iPlan_) {
         cufftDestroy(iPlan_);
      }
   }

   /*
   * Setup mesh dimensions.
   */
   template <int D>
   void FFT<D>::setup(IntVec<D> const& meshDimensions)
   {
      // Precondition
      UTIL_CHECK(!isSetup_);

      // Create local r-grid and k-grid field objects
      RField<D> rField;
      rField.allocate(meshDimensions);
      RFieldDft<D> kField;
      kField.allocate(meshDimensions);

      setup(rField, kField);
   }

   /*
   * Check and (if necessary) setup mesh dimensions.
   */
   template <int D>
   void FFT<D>::setup(RField<D>& rField, RFieldDft<D>& kField)
   {
      // Preconditions
      UTIL_CHECK(!isSetup_);
      IntVec<D> rDimensions = rField.meshDimensions();
      IntVec<D> kDimensions = kField.meshDimensions();
      UTIL_CHECK(rDimensions == kDimensions);

      // Set Mesh dimensions
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(rDimensions[i] > 0);
         meshDimensions_[i] = rDimensions[i];
         rSize_ *= rDimensions[i];
         if (i < D - 1) {
            kSize_ *= rDimensions[i];
         } else {
            kSize_ *= (rDimensions[i]/2 + 1);
         }
      }

      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);


      // Make FFTW plans (explicit specializations)
      makePlans(rField, kField);

      // Allocate rFieldCopy_ array if necessary
      if (!rFieldCopy_.isAllocated()) {
          rFieldCopy_.allocate(rDimensions);
      } else {
          if (rFieldCopy_.capacity() != rSize_) {
             rFieldCopy_.deallocate();
             rFieldCopy_.allocate(rDimensions);
          }
      }
      UTIL_CHECK(rFieldCopy_.capacity() == rSize_);

      // Allocate kFieldCopy_ array if necessary
      if (!kFieldCopy_.isAllocated()) {
          kFieldCopy_.allocate(kDimensions);
      } else {
          if (kFieldCopy_.capacity() != rSize_) {
             kFieldCopy_.deallocate();
             kFieldCopy_.allocate(kDimensions);
          }
      }
      UTIL_CHECK(kFieldCopy_.capacity() == kSize_);

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
      scaleRealData<<<nBlocks, nThreads>>>(rField.cField(), scale, rSize_);
      
      //perform fft
      #ifdef SINGLE_PRECISION
      if(cufftExecR2C(fPlan_, rField.cField(), kField.cField()) != CUFFT_SUCCESS) {
         std::cout << "CUFFT error: forward" << std::endl;
         return;
      }
      #else
      if(cufftExecD2Z(fPlan_, rField.cField(), kField.cField()) != CUFFT_SUCCESS) {
         std::cout << "CUFFT error: forward" << std::endl;
         return;
      }
      #endif

   }

   /*
   * Execute forward transform without destroying input.
   */
   template <int D>
   void FFT<D>::forwardTransformSafe(RField<D> const & rField, RFieldDft<D>& kField)
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

      #ifdef SINGLE_PRECISION
      if(cufftExecC2R(iPlan_, kField.cField(), rField.cField()) != CUFFT_SUCCESS) {
         std::cout << "CUFFT error: inverse" << std::endl;
         return;
      }
      #else
      if(cufftExecZ2D(iPlan_, kField.cField(), rField.cField()) != CUFFT_SUCCESS) {
         std::cout << "CUFFT error: inverse" << std::endl;
         return;
      }
      #endif
   
   }

   /*
   * Execute inverse (complex-to-real) transform without destroying input.
   */
   template <int D>
   void FFT<D>::inverseTransformSafe(RFieldDft<D> const & kField, RField<D>& rField) 
   const
   {
      UTIL_CHECK(kFieldCopy_.capacity()==kField.capacity());

      kFieldCopy_ = kField;
      inverseTransform(kFieldCopy_, rField);
   }

}
}

#if 0
static __global__ void scaleRealData(cudaReal* data, cudaReal scale, int size) {
   
   //write code that will scale
   int nThreads = blockDim.x * gridDim.x;
   int startId = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startId; i < size; i += nThreads ) {
      data[i] *= scale;
   }
   
}
#endif

#endif
