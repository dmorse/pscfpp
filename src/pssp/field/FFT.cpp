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

   // Forward transform, explicit specializations.

   template<>
   void FFT<1>::forwardTransform(RField<1>& in, RFieldDFT<1>& out)
   {
      setup(in, out);
      unsigned int flags = FFTW_ESTIMATE;
      if (fPlan_ == 0) {
         fPlan_ = fftw_plan_dft_r2c_1d(rSize_, &in[0], &out[0], flags);
         fftw_execute(fPlan_);
      }
   }

   template<>
   void FFT<2>::forwardTransform(RField<2>& in, RFieldDFT<2>& out)
   {
      setup(in, out);
   }

   template<>
   void FFT<3>::forwardTransform(RField<3>& in, RFieldDFT<3>& out)
   {
      setup(in, out);
   }

}
