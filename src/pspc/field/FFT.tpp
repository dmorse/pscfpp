#ifndef PSPC_FFT_TPP
#define PSPC_FFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

namespace Pscf {
namespace Pspc
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
      UTIL_CHECK(rDimensions == kField.meshDimensions());

      // Set and check mesh dimensions
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

      // Allocate work array if necessary
      if (!work_.isAllocated()) {
          work_.allocate(rDimensions);
      } else {
          if (work_.capacity() != rSize_) {
             work_.deallocate();
             work_.allocate(rDimensions);
          }
      }
      UTIL_CHECK(work_.capacity() == rSize_);

      // Make FFTW plans (see explicit specializations FFT.cpp)
      makePlans(rField, kField);

      isSetup_ = true;
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(RField<D> const & rField, RFieldDft<D>& kField)   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(work_.capacity() == rSize_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Copy rescaled input data prior to work array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < rSize_; ++i) {
         work_[i] = rField[i]*scale;
      }
     
      // Execute preplanned forward transform 
      fftw_execute_dft_r2c(fPlan_, &work_[0], &kField[0]);
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void 
   FFT<D>::inverseTransform(RFieldDft<D> & kField, RField<D>& rField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      fftw_execute_dft_c2r(iPlan_, &kField[0], &rField[0]);

   }

}
}
#endif
