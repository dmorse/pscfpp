#ifndef PRDC_CPU_FFT_TPP
#define PRDC_CPU_FFT_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFT<D>::FFT()
    : rFieldCopy_(),
      meshDimensions_(0),
      rSize_(0),
      kSize_(0),
      rcfPlan_(0),
      criPlan_(0),
      ccfPlan_(0),
      cciPlan_(0),
      isSetup_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FFT<D>::~FFT()
   {
      if (rcfPlan_) {
         fftw_destroy_plan(rcfPlan_);
      }
      if (criPlan_) {
         fftw_destroy_plan(criPlan_);
      }
      if (ccfPlan_) {
         fftw_destroy_plan(ccfPlan_);
      }
      if (cciPlan_) {
         fftw_destroy_plan(cciPlan_);
      }
   }

   /*
   * Setup mesh dimensions, work memory and FFT plans.
   */
   template <int D>
   void FFT<D>::setup(IntVec<D> const& meshDimensions)
   {
      // Precondition
      UTIL_CHECK(!isSetup_);

      // Set and check mesh dimensions
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
      UTIL_CHECK(meshDimensions == rFieldCopy_.meshDimensions());
      UTIL_CHECK(rFieldCopy_.capacity() == rSize_);

      // Allocate kFieldCopy_ array if necessary
      if (!kFieldCopy_.isAllocated()) {
          kFieldCopy_.allocate(meshDimensions);
      } else {
          if (kFieldCopy_.capacity() != kSize_) {
             kFieldCopy_.deallocate();
             kFieldCopy_.allocate(meshDimensions);
          }
      }
      UTIL_CHECK(meshDimensions == kFieldCopy_.meshDimensions());
      UTIL_CHECK(kFieldCopy_.capacity() == kSize_);

      // Allocate cFieldCopy_ array if necessary
      if (!cFieldCopy_.isAllocated()) {
          cFieldCopy_.allocate(meshDimensions);
      } else {
          if (cFieldCopy_.capacity() != rSize_) {
             cFieldCopy_.deallocate();
             cFieldCopy_.allocate(meshDimensions);
          }
      }
      UTIL_CHECK(meshDimensions == cFieldCopy_.meshDimensions());
      UTIL_CHECK(cFieldCopy_.capacity() == rSize_);

      #if 0
      // Create local field objects used for plans
      RField<D> rField;
      rField.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == rField.meshDimensions());
      UTIL_CHECK(rField.capacity() == rSize_);

      RFieldDft<D> kField;
      kField.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == kField.meshDimensions());
      UTIL_CHECK(kField.capacity() == kSize_);

      CField<D> cFieldIn;
      cFieldIn.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == cFieldIn.meshDimensions());
      UTIL_CHECK(cFieldIn.capacity() == rSize_);
      #endif

      CField<D> cFieldOut;
      cFieldOut.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == cFieldOut.meshDimensions());
      UTIL_CHECK(cFieldOut.capacity() == rSize_);

      // Make FFTW plans (see explicit specializations FFT.cpp)
      //makePlans(rField, kField, cFieldIn, cFieldOut);
      makePlans(rFieldCopy_, kFieldCopy_, cFieldCopy_, cFieldOut);

      isSetup_ = true;
   }

   // Real <-> Complex Transforms

   /*
   * Execute real-to-complex forward transform.
   */
   template <int D>
   void FFT<D>::forwardTransform(RField<D> const & rField, 
                                 RFieldDft<D>& kField)   
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rFieldCopy_.capacity() == rSize_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Copy rescaled input data prior to work array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < rSize_; ++i) {
         rFieldCopy_[i] = rField[i]*scale;
      }
     
      // Execute preplanned forward transform 
      fftw_execute_dft_r2c(rcfPlan_, &rFieldCopy_[0], &kField[0]);
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

      fftw_execute_dft_c2r(criPlan_, &kField[0], &rField[0]);

   }

   /*
   * Execute inverse (complex-to-real) transform without destroying input.
   */
   template <int D>
   void 
   FFT<D>::inverseTransformSafe(RFieldDft<D> const & kField, 
                                RField<D>& rField) 
   const
   {
      UTIL_CHECK(kFieldCopy_.capacity()==kField.capacity());

      kFieldCopy_ = kField;
      inverseTransform(kFieldCopy_, rField);
   }

   // Complex <-> Complex Transforms

   /*
   * Execute complex-to-complex forward transform.
   */
   template <int D>
   void 
   FFT<D>::forwardTransform(CField<D> const & rField, CField<D>& kField)   
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(cFieldCopy_.capacity() == rSize_);
      UTIL_CHECK(cFieldCopy_.meshDimensions() == meshDimensions_);
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Copy rescaled input data prior to work array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < rSize_; ++i) {
         cFieldCopy_[i][0] = rField[i][0]*scale;
         cFieldCopy_[i][1] = rField[i][1]*scale;
      }
     
      // Execute preplanned forward transform 
      fftw_execute_dft(ccfPlan_, &cFieldCopy_[0], &kField[0]);
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void 
   FFT<D>::inverseTransform(CField<D> & kField, CField<D>& rField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      fftw_execute_dft(cciPlan_, &kField[0], &rField[0]);
   }

}
}
}
#endif
