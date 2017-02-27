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
      rPlan_(0),
      hasDimensions_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FFT<D>::~FFT()
   {
      if (fPlan_) {
         fPlan_ = 0;
      }
      if (rPlan_) {
         rPlan_ = 0;
      }
   }

   /*
   * Check and (if necessary) setup mesh dimensions.
   */
   template <int D>
   void FFT<D>::setup(const RField<D>& rField, const RFieldDFT<D>& kField)
   {
      IntVec<D> rDimensions = rField.meshDimensions();
      IntVec<D> kDimensions = kField.meshDimensions();
      UTIL_CHECK(rDimensions == kDimensions);
      if (hasDimensions_) {
         UTIL_CHECK(rDimensions == meshDimensions_);
      } else {
         setDimensions(rDimensions);
      }
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
   }

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void FFT<D>::setDimensions(const IntVec<D>& meshDimensions)
   {
      UTIL_CHECK(!hasDimensions_);
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
      hasDimensions_ = true;
   }

}
#endif
