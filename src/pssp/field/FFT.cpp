/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

namespace Pscf {
namespace Pssp {

   using namespace Util;

   // Forward transform, explicit specializations.

   template<>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField)
   {
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_1d(rSize_, &rField[0], &kField[0], flags);
      iPlan_ = fftw_plan_dft_c2r_1d(rSize_, &kField[0], &rField[0], flags);
   }

}
}
