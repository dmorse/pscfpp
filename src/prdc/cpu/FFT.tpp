#ifndef PRDC_CPU_FFT_TPP
#define PRDC_CPU_FFT_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

/*
* A note about const_casts:
* 
* The FFTW library is used in this file to perform discrete Fourier 
* transforms. FFTW's complex-to-real inverse transform overwrites its
* input array, but all other out-of-place transforms leave the input 
* array unaltered. However, all transforms in the FFTW library require
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
* with the FFTW library's requirement of non-const inputs. This conflict 
* is resolved using a const_cast, in which the const pointer to the input 
* array is made non-const when passed into FFTW functions. The use of
* const_cast is reserved only for these few select cases in which we are
* confident that the input array will not be modified.
*
* For more information about the relevant FFTW methods, see the FFTW 
* documentation at https://www.fftw.org/fftw3_doc/index.html. In Section 
* 4.3.2, it is stated that, by default, "an out-of-place transform must 
* not change its input array," except for complex-to-real transforms, in 
* which case "no input-preserving algorithms are implemented." Finally, 
* we note that the unit tests for this FFT class check that the input 
* array is unaltered, allowing developers to continually ensure that 
* the FFTW functions do not modify their input unexpectedly.
*/

namespace Pscf {
namespace Prdc {
namespace Cpu {

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
   void FFT<D>::setup(IntVec<D> const & meshDimensions)
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

      // Create temporary RField and CField objects used for plans
      RField<D> rField;
      rField.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == rField.meshDimensions());
      UTIL_CHECK(rField.capacity() == rSize_);

      CField<D> cFieldIn;
      cFieldIn.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == cFieldIn.meshDimensions());
      UTIL_CHECK(cFieldIn.capacity() == rSize_);

      CField<D> cFieldOut;
      cFieldOut.allocate(meshDimensions);
      UTIL_CHECK(meshDimensions == cFieldOut.meshDimensions());
      UTIL_CHECK(cFieldOut.capacity() == rSize_);

      // Make FFTW plans (see explicit specializations in FFT.cpp)
      makePlans(rField, kFieldCopy_, cFieldIn, cFieldOut);

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
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);
     
      // Execute preplanned forward transform 
      // (See note at top of file explaining this use of const_cast)
      fftw_execute_dft_r2c(rcfPlan_, const_cast<double*>(&rField[0]), 
                           &kField[0]);

      // Rescale the resulting array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < kSize_; ++i) {
         kField[i][0] *= scale;
         kField[i][1] *= scale;
      }
   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void 
   FFT<D>::inverseTransformUnsafe(RFieldDft<D> & kField, RField<D>& rField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Execute preplanned inverse transform 
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
      UTIL_CHECK(kFieldCopy_.capacity() == kField.capacity());

      kFieldCopy_ = kField;
      inverseTransformUnsafe(kFieldCopy_, rField);
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
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);
     
      // Execute preplanned forward transform
      // (See note at top of file explaining this use of const_cast) 
      fftw_execute_dft(ccfPlan_, const_cast<fftw_complex*>(&rField[0]), 
                       &kField[0]);

      // Rescale the resulting array
      double scale = 1.0/double(rSize_);
      for (int i = 0; i < rSize_; ++i) {
         kField[i][0] *= scale;
         kField[i][1] *= scale;
      }
   }

   /*
   * Execute inverse (complex-to-complex) transform.
   */
   template <int D>
   void 
   FFT<D>::inverseTransform(CField<D> const & kField, CField<D>& rField)
   const
   {
      UTIL_CHECK(isSetup_)
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == rSize_);
      UTIL_CHECK(rField.meshDimensions() == meshDimensions_);
      UTIL_CHECK(kField.meshDimensions() == meshDimensions_);

      // Execute preplanned inverse transform
      // (See note at top of file explaining this use of const_cast) 
      fftw_execute_dft(cciPlan_, const_cast<fftw_complex*>(&kField[0]), 
                       &rField[0]);
   }

   // Static function

   /*
   * Compute dimensions and size of the k-size mesh (static)
   */
   template <int D>
   void FFT<D>::computeKMesh(IntVec<D> const & rMeshDimensions,
                             IntVec<D> & kMeshDimensions,
                             int & kSize )
   {
      kSize = 1;
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions[i] = rMeshDimensions[i];
         } else {
            kMeshDimensions[i] = rMeshDimensions[i]/2 + 1;
         }
         kSize *= kMeshDimensions[i];
      }
   }

}
}
}
#endif
