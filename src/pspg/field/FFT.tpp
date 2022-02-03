#ifndef PSPG_FFT_TPP
#define PSPG_FFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"
#include <pspg/GpuResources.h>

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
   * Check and (if necessary) setup mesh dimensions.
   */
   template <int D>
   void FFT<D>::setup(RDField<D>& rField, RDFieldDft<D>& kField)
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

      isSetup_ = true;
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(RDField<D>& rField, RDFieldDft<D>& kField)
   {
      // Check dimensions or setup
      if (isSetup_) {
         UTIL_CHECK(rField.capacity() == rSize_);
         UTIL_CHECK(kField.capacity() == kSize_);
      } else {
         UTIL_THROW("Reached the spot it should 'never reach'.");
         //should never reach here in a parallel block. breaks instantly
         setup(rField, kField);
      }

      // Copy rescaled input data prior to work array
      cudaReal scale = 1.0/cudaReal(rSize_);
      scaleRealData<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(rField.cDField(), scale, rSize_);
      
      //perform fft
      #ifdef SINGLE_PRECISION
      if(cufftExecR2C(fPlan_, rField.cDField(), kField.cDField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: forward"<<std::endl;
         return;
      }
      #else
      if(cufftExecD2Z(fPlan_, rField.cDField(), kField.cDField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: forward"<<std::endl;
         return;
      }
      #endif

   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFT<D>::inverseTransform(RDFieldDft<D>& kField, RDField<D>& rField)
   {
      if (!isSetup_) {
         //UTIL_CHECK(0);
         setup(rField, kField);
      }

      #ifdef SINGLE_PRECISION
      if(cufftExecC2R(iPlan_, kField.cDField(), rField.cDField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: inverse"<<std::endl;
         return;
      }
      #else
      if(cufftExecZ2D(iPlan_, kField.cDField(), rField.cDField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: inverse"<<std::endl;
         return;
      }
      #endif
   
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
