#ifndef FFT_TPP
#define FFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

namespace Pscf {
namespace Pssp
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFT<D>::FFT()
    : work_(),
      meshDimensions_(0),
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
      work_.allocate(rDimensions);
      UTIL_CHECK(work_.capacity() == rSize_);
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
   void FFT<D>::forwardTransform(RField<D>& rField, RFieldDft<D>& kField)
   {
      // Check dimensions or setup
      if (isSetup_) {
         UTIL_CHECK(work_.capacity() == rSize_);
         UTIL_CHECK(rField.capacity() == rSize_);
         UTIL_CHECK(kField.capacity() == kSize_);
      } else {
         setup(rField, kField);
      }

      // Copy rescaled input data prior to work array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < rSize_; ++i) {
         work_[i] = rField[i]*scale;
      }
      
      //Are there any instances where the original array is important?
      fftw_execute_dft_r2c(fPlan_, &work_[0], &kField[0]);
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFT<D>::inverseTransform(RFieldDft<D>& kField, RField<D>& rField)
   {
      if (!isSetup_) {
         setup(rField, kField);
         fftw_execute(iPlan_);
      } else {
         fftw_execute_dft_c2r(iPlan_, &kField[0], &rField[0]);
      }
   }

}
}
#endif
