#ifndef PRDC_CUDA_FFT_BATCHED_TPP
#define PRDC_CUDA_FFT_BATCHED_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.h"
#include <pscf/cuda/GpuResources.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFTBatched<D>::FFTBatched()
    : meshDimensions_(0),
      kMeshDimensions_(0),
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
   FFTBatched<D>::~FFTBatched()
   {
      if (fPlan_) {
         cufftDestroy(fPlan_);
      }
      if (iPlan_) {
         cufftDestroy(iPlan_);
      }
   }

   /*
   * Set up FFT calculation (store grid dimensions and make FFT plan)
   */
   template <int D>
   void FFTBatched<D>::setup(const IntVec<D>& rDim, const IntVec<D>& kDim, 
                             int batchSize)
   {
      // Preconditions
      UTIL_CHECK(!isSetup_);

      // Set Mesh dimensions
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(rDim[i] > 0);
         meshDimensions_[i] = rDim[i];
         kMeshDimensions_[i] = kDim[i];
         rSize_ *= rDim[i];
         if (i < D - 1) {
            kSize_ *= rDim[i];
         } else {
            kSize_ *= (rDim[i]/2 + 1);
         }
      }

      // Make FFT plans
      makePlans(rDim, kDim, batchSize);

      isSetup_ = true;
   }

   /**
   * Set the batch size to a new value. isSetup() must already be true.
   */
   template <int D>
   void FFTBatched<D>::setBatchSize(int batchSize)
   {
      UTIL_CHECK(isSetup_);

      if (batchSize == batchSize_) {
         // Nothing to do
         return; 
      } else {
         // Remake FFT plans
         makePlans(meshDimensions_, kMeshDimensions_, batchSize); 
      }
   }

   /*
   * Make plans for variable batch size
   */
   template <int D>
   void FFTBatched<D>::makePlans(const IntVec<D>& rDim, 
                                 const IntVec<D>& kDim, int batchSize)
   {
      batchSize_ = batchSize;

      int n[D];
      for(int i = 0; i < D; i++) {
         n[i] = rDim[i];
      }

      #ifdef SINGLE_PRECISION
      if (cufftPlanMany(&fPlan_, D, n,   //plan, rank, n
                        NULL, 1, rSize_, //inembed, istride, idist
                        NULL, 1, kSize_, //onembed, ostride, odist
                        CUFFT_R2C, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if (cufftPlanMany(&iPlan_, D, n,   //plan, rank, n
                        NULL, 1, kSize_, //inembed, istride, idist
                        NULL, 1, rSize_, //onembed, ostride, odist
                        CUFFT_C2R, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      #else
      if (cufftPlanMany(&fPlan_, D, n,   // plan, rank, n
                        NULL, 1, rSize_, // inembed, istride, idist
                        NULL, 1, kSize_, // onembed, ostride, odist
                        CUFFT_D2Z, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if (cufftPlanMany(&iPlan_, D, n,   // plan, rank, n
                        NULL, 1, kSize_, // inembed, istride, idist
                        NULL, 1, rSize_, // onembed, ostride, odist
                        CUFFT_Z2D, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      #endif
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFTBatched<D>::forwardTransform(DeviceArray<cudaReal>& in, 
                                        DeviceArray<cudaComplex>& out) const
   {
      int nr = rSize_ * batchSize_;
      int nk = kSize_ * batchSize_;
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(in.capacity() == nr);
      UTIL_CHECK(out.capacity() == nk);

      // Scale for every batch
      cudaReal scale = 1.0/cudaReal(rSize_);
      VecOp::mulEqS(in, scale);
      
      // Perform FFT
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecR2C(fPlan_, in.cArray(), out.cArray());
      #else
      result = cufftExecD2Z(fPlan_, in.cArray(), out.cArray());
      #endif

      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft real-to-complex forward transform.");
      }
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFTBatched<D>::inverseTransform(DeviceArray<cudaComplex>& in, 
                                        DeviceArray<cudaReal>& out) const
   {
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(in.capacity() == kSize_ * batchSize_);
      UTIL_CHECK(out.capacity() == rSize_ * batchSize_);

      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2R(iPlan_, in.cArray(), out.cArray());
      #else
      result = cufftExecZ2D(iPlan_, in.cArray(), out.cArray());
      #endif

      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft complex-to-real inverse transform.");
      }
   }

}
}
}
#endif
