#ifndef PRDC_CUDA_FFT_BATCHED_TPP
#define PRDC_CUDA_FFT_BATCHED_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.h"
#include "VecOp.h"

/*
* A note about const_casts:
* 
* The cuFFT library is used in this file to perform discrete Fourier 
* transforms. cuFFT's complex-to-real inverse transform overwrites its
* input array, but all other out-of-place transforms leave the input 
* array unaltered. However, all transforms in the cuFFT library require
* non-const pointers to the input array, even though they do not alter 
* the array.
*
* In order to maintain const-correctness in PSCF, this FFT class accepts
* const input arrays for its methods that perform a Fourier transform,
* unless the transform is expected to modify / overwrite its input (as
* is the case for complex-to-real inverse transforms). This denotes to
* the caller of the method that the input array will not be altered,
* which is an accurate representation of the expected behavior.
* 
* However, the const-correctness of this FFT class creates a conflict
* with the cuFFT library's requirement of non-const inputs. This conflict 
* is resolved using a const_cast, in which the const pointer to the input 
* array is made non-const when passed into cuFFT functions. The use of
* const_cast is reserved only for these few select cases in which we are
* confident that the input array will not be modified.
*
* For more information about the relevant cuFFT methods, see the cuFFT 
* documentation at https://docs.nvidia.com/cuda/cufft/. Unfortunately,
* there is no explicitly documented guarantee that the transforms do not
* modify their input, though it is implied. In Section 2.4, it is stated
* that "[o]ut-of-place complex-to-real FFT will always overwrite input 
* buffer." No such claim is made for any other out-of-place transforms, 
* implying that they do not overwrite their inputs. Further, "[t]he 
* cuFFT API is modeled after FFTW," (beginning of Section 2), and FFTW
* is much more explicit in their documentation 
* (https://www.fftw.org/fftw3_doc/index.html). In Section 4.3.2, it is
* stated that, by default, "an out-of-place transform must not change
* its input array," except for complex-to-real transforms, in which case
* "no input-preserving algorithms are implemented." Finally, we note 
* that the unit tests for this FFT class check that the input array is
* unaltered, allowing developers to continually ensure that the cuFFT
* functions do not modify their input unexpectedly.
*/

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
   void FFTBatched<D>::setup(const IntVec<D>& meshDimensions, int batchSize)
   {
      // Preconditions
      UTIL_CHECK(!isSetup_);

      // Set Mesh dimensions
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         if (i < D - 1) {
            kMeshDimensions_[i] = meshDimensions[i];
         } else {
            kMeshDimensions_[i] = (meshDimensions[i]/2 + 1);
         }
         rSize_ *= meshDimensions_[i];
         kSize_ *= kMeshDimensions_[i];
      }

      // Make FFT plans
      makePlans(batchSize);

      isSetup_ = true;
   }

   /**
   * Set the batch size to a new value. isSetup() must already be true.
   */
   template <int D>
   void FFTBatched<D>::resetBatchSize(int batchSize)
   {
      UTIL_CHECK(isSetup_);

      if (batchSize == batchSize_) {
         // Nothing to do
         return; 
      } else {
         // Remake FFT plans
         makePlans(batchSize); 
      }
   }

   /*
   * Make plans for variable batch size
   */
   template <int D>
   void FFTBatched<D>::makePlans(int batchSize)
   {
      batchSize_ = batchSize;

      UTIL_CHECK(kSize_ > 0);
      UTIL_CHECK(rSize_ > 0);

      int n[D];
      for(int i = 0; i < D; i++) {
         UTIL_CHECK(meshDimensions_[i] > 0);
         n[i] = meshDimensions_[i];
      }

      #ifdef SINGLE_PRECISION
      if (cufftPlanMany(&fPlan_, D, n,   //plan, rank, n
                        NULL, 1, rSize_, //inembed, istride, idist
                        NULL, 1, kSize_, //onembed, ostride, odist
                        CUFFT_R2C, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"FFTBatched: plan creation failed "<<std::endl;
         exit(1);
      }
      if (cufftPlanMany(&iPlan_, D, n,   //plan, rank, n
                        NULL, 1, kSize_, //inembed, istride, idist
                        NULL, 1, rSize_, //onembed, ostride, odist
                        CUFFT_C2R, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"FFTBatched: plan creation failed "<<std::endl;
         exit(1);
      }
      #else
      if (cufftPlanMany(&fPlan_, D, n,   // plan, rank, n
                        NULL, 1, rSize_, // inembed, istride, idist
                        NULL, 1, kSize_, // onembed, ostride, odist
                        CUFFT_D2Z, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"FFTBatched: plan creation failed "<<std::endl;
         exit(1);
      }
      if (cufftPlanMany(&iPlan_, D, n,   // plan, rank, n
                        NULL, 1, kSize_, // inembed, istride, idist
                        NULL, 1, rSize_, // onembed, ostride, odist
                        CUFFT_Z2D, batchSize_) != CUFFT_SUCCESS) {
         std::cout<<"FFTBatched: plan creation failed "<<std::endl;
         exit(1);
      }
      #endif
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFTBatched<D>::forwardTransform(DeviceArray<cudaReal> const & rFields, 
                                        DeviceArray<cudaComplex>& kFields) 
   const
   {
      // Preconditions
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rFields.capacity() == rSize_ * batchSize_);
      UTIL_CHECK(kFields.capacity() == kSize_ * batchSize_);
      
      // Perform FFT
      // (See note at top of file explaining this use of const_cast)
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecR2C(fPlan_, const_cast<cudaReal*>(rFields.cArray()), 
                            kFields.cArray());
      #else
      result = cufftExecD2Z(fPlan_, const_cast<cudaReal*>(rFields.cArray()), 
                            kFields.cArray());
      #endif

      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft real-to-complex forward transform.");
      }

      // Rescale output data in-place
      cudaReal scale = 1.0/cudaReal(rSize_);
      VecOp::mulEqS(kFields, scale);
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void 
   FFTBatched<D>::inverseTransformUnsafe(DeviceArray<cudaComplex>& kFields, 
                                         DeviceArray<cudaReal>& rFields) 
   const
   {
      // Preconditions
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(kFields.capacity() == kSize_ * batchSize_);
      UTIL_CHECK(rFields.capacity() == rSize_ * batchSize_);

      // Perform FFT
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2R(iPlan_, kFields.cArray(), rFields.cArray());
      #else
      result = cufftExecZ2D(iPlan_, kFields.cArray(), rFields.cArray());
      #endif

      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft complex-to-real inverse transform.");
      }
   }

}
}
}
#endif
