#ifndef FFT_TPP
#define FFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

namespace Pssp
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
         fftw_destroy_plan(fPlan_);
      }
      if (iPlan_) {
         fftw_destroy_plan(iPlan_);
      }
   }

   /*
   * Check and (if necessary) setup mesh dimensions.
   */
   template <int D>
   void FFT<D>::setup(RField<D>& rField, RFieldDFT<D>& kField)
   {
      IntVec<D> rDimensions = rField.meshDimensions();
      IntVec<D> kDimensions = kField.meshDimensions();
      UTIL_CHECK(rDimensions == kDimensions);
      if (isSetup_) {
         UTIL_CHECK(rDimensions == meshDimensions_);
         UTIL_CHECK(rField.capacity() == rSize_);
         UTIL_CHECK(kField.capacity() == kSize_);
      } else {
         setDimensions(rDimensions);
         UTIL_CHECK(rField.capacity() == rSize_);
         UTIL_CHECK(kField.capacity() == kSize_);
         makePlans(rField, kField);
         isSetup_ = true;
      }
   }

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void FFT<D>::setDimensions(const IntVec<D>& meshDimensions)
   {
      int rSize_ = 1;
      int kSize_ = 1;
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
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(RField<D>& rField, RFieldDFT<D>& kField)
   {
      if (!isSetup_) {
         setup(rField, kField);
         fftw_execute(fPlan_);
      } else {
         fftw_execute_dft_r2c(fPlan_, &rField[0], &kField[0]);
      }
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFT<D>::inverseTransform(RFieldDFT<D>& kField, RField<D>& rField)
   {
      if (!isSetup_) {
         setup(rField, kField);
         fftw_execute(iPlan_);
      } else {
         fftw_execute_dft_c2r(fPlan_, &kField[0], &rField[0]);
      }
   }

}
#endif
