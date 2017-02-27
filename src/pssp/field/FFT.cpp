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
   FFT::FFT()
    : meshDimensions_(0),
      spaceDimension_(0),
      rSize_(0),
      kSize_(0),
      fPlan_(0),
      rPlan_(0)
   {}

   /*
   * Destructor.
   */
   FFT::~FFT()
   {
      if (fPlan_) {
         fPlan_ = 0;
      }
      if (rPlan_) {
         rPlan_ = 0;
      }
   }

   // Forward transform, explicit specializations.

   template<>
   void FFT::forwardTransform<1>(RMeshField& in, KMeshField& out)
   {
      setup<1>(in, out);
      unsigned int flags = FFTW_ESTIMATE;
      if (fPlan_ == 0) {
         fPlan_ = fftw_plan_dft_r2c_1d(rSize_, &in[0], &out[0], flags);
         fftw_execute(fPlan_);
      }
   }

   template<>
   void FFT::forwardTransform<2>(RMeshField& in, KMeshField& out)
   {
      setup<2>(in, out);
   }

   template<>
   void FFT::forwardTransform<3>(RMeshField& in, KMeshField& out)
   {
      setup<3>(in, out);
   }

}
