/*
* PSCF++ Package
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   // Explicit class instantiations

   template class FFT<1>;
   template class FFT<2>;
   template class FFT<3>;


   // Planning functions, explicit specializations.

   template<>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField)
   {
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_1d(rSize_, &rField[0], &kField[0], flags);
      iPlan_ = fftw_plan_dft_c2r_1d(rSize_, &kField[0], &rField[0], flags);
   }

   template <>
   void FFT<2>::makePlans(RField<2>& rField, RFieldDft<2>& kField)
   {
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_2d(meshDimensions_[0], meshDimensions_[1],
      	                           &rField[0], &kField[0], flags);
      iPlan_ = fftw_plan_dft_c2r_2d(meshDimensions_[0], meshDimensions_[1],
                                    &kField[0], &rField[0], flags);
   }

   template <>
   void FFT<3>::makePlans(RField<3>& rField, RFieldDft<3>& kField)
   {
      unsigned int flags = FFTW_ESTIMATE;
      fPlan_ = fftw_plan_dft_r2c_3d(meshDimensions_[0], meshDimensions_[1],
      	                           meshDimensions_[2], &rField[0], &kField[0],
      	                           flags);
      iPlan_ = fftw_plan_dft_c2r_3d(meshDimensions_[0], meshDimensions_[1],
                                    meshDimensions_[2], &kField[0], &rField[0],
                                    flags);
   }

}
}
