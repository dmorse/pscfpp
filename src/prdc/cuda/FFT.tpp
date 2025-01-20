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
   * Setup grid dimensions, plans and work space.
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

      // Reallocate kFieldCopy_ array if necessary
      if (kFieldCopy_.isAllocated()) {
         if (kFieldCopy_.capacity() != kSize_) {
            kFieldCopy_.deallocate();
            kFieldCopy_.allocate(meshDimensions);
         }
      }

      // Make FFTW plans (explicit specializations)
      makePlans();

      isSetup_ = true;
   }

   /*
   * Compute forward (real-to-complex) discrete Fourier transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(RField<D> const & rField, 
                                 RFieldDft<D>& kField) const
   {
      // Preconditions
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      
      // Perform transform
      // (See note at top of file explaining this use of const_cast)
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecR2C(rcfPlan_, const_cast<cudaReal*>(rField.cArray()), 
                            kField.cArray());
      #else
      result = cufftExecD2Z(rcfPlan_, const_cast<cudaReal*>(rField.cArray()), 
                            kField.cArray());
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft real-to-complex forward transform");
      }

      // Rescale output data in-place
      cudaReal scale = 1.0/cudaReal(rSize_);
      VecOp::mulEqS(kField, scale);
   }

   /*
   * Compute inverse (complex-to-real) DFT, overwriting the input.
   */
   template <int D>
   void FFT<D>::inverseTransformUnsafe(RFieldDft<D>& kField, 
                                       RField<D>& rField) const
   {
      // Preconditions
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
   * Compute inverse (complex-to-real) DFT without overwriting input.
   */
   template <int D>
   void FFT<D>::inverseTransformSafe(RFieldDft<D> const & kField, 
                                     RField<D>& rField) const
   {
      // if kFieldCopy_ has been previously allocated, check size is correct
      if (kFieldCopy_.isAllocated()) { 
         UTIL_CHECK(kFieldCopy_.capacity() == kSize_);
         UTIL_CHECK(kFieldCopy_.meshDimensions() == meshDimensions_);
      }

      // make copy of kField (allocates kFieldCopy_ if necessary)
      kFieldCopy_ = kField; 

      // Perform transform using copy of kField
      inverseTransformUnsafe(kFieldCopy_, rField);
   }

   // Complex-to-Complex Transforms

   /*
   * Execute forward complex-to-complex transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(CField<D> const & rField, CField<D>& kField)
   const
   {
      // Preconditions
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);
      
      // Perform transform
      // (See note at top of file explaining this use of const_cast)
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2C(ccPlan_, const_cast<cudaComplex*>(rField.cArray()), 
                            kField.cArray(), CUFFT_FORWARD);
      #else
      result = cufftExecZ2Z(ccPlan_, const_cast<cudaComplex*>(rField.cArray()), 
                            kField.cArray(), CUFFT_FORWARD);
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW("Failure in cufft complex-to-complex forward transform");
      }

      // Rescale output data in-place
      cudaReal scale = 1.0/cudaReal(rSize_);
      VecOp::mulEqS(kField, scale);
   }

   /*
   * Execute inverse (complex-to-complex) transform.
   */
   template <int D>
   void FFT<D>::inverseTransform(CField<D> const & kField, CField<D>& rField) 
   const
   {
      // Preconditions
      UTIL_CHECK(isSetup_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Perform transform
      // (See note at top of file explaining this use of const_cast)
      cufftResult result;
      #ifdef SINGLE_PRECISION
      result = cufftExecC2C(ccPlan_, const_cast<cudaComplex*>(kField.cArray()), 
                            rField.cArray(), CUFFT_INVERSE);
      #else
      result = cufftExecZ2Z(ccPlan_, const_cast<cudaComplex*>(kField.cArray()), 
                            rField.cArray(), CUFFT_INVERSE);
      #endif
      if (result != CUFFT_SUCCESS) {
         UTIL_THROW( "Failure in cufft complex-to-complex inverse transform");
      }
   
   }

}
}
}
#endif
