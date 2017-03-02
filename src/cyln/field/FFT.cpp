/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
   FFT::FFT()
    : work_(),
      rSize_(0),
      kSize_(0),
      fPlan_(0),
      iPlan_(0),
      isSetup_(false)
   {}

   /*
   * Destructor.
   */
   FFT::~FFT()
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
   void FFT::setup(Array<double>& rField, Array<fftw_complex>& kField)
   {
      // Preconditions
      UTIL_CHECK(!isSetup_);
      IntVec<D> rSize_ = rField.capacity();
      IntVec<D> kSize_ = kField.capacity();
      UTIL_CHECK(kSize_ == rSize_/2 + 1);

      work_.allocate(rSize_);

      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_1d(rSize_, &rField[0], &kField[0], flags);
      iPlan_ = fftw_plan_dft_c2r_1d(rSize_, &kField[0], &rField[0], flags);

      isSetup_ = true;
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFT::forwardTransform(Array<double>& rField, 
                              Array<fftw_complex>& kField)
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

      fftw_execute_dft_r2c(fPlan_, &work_[0], &kField[0]);
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   void FFT::inverseTransform(Array<fftw_complex>& kField, 
                              Array<double>& rField)
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
