/*
* PSCF Package
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   // Explicit class instantiations

   template class FFT<1>;
   template class FFT<2>;
   template class FFT<3>;

   // Planning functions, explicit specializations.

   template<>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField, 
                          CField<1>& cFieldIn, CField<1>& cFieldOut)
   {
      int n0 = rSize_;
      unsigned int flags = FFTW_ESTIMATE;
      rcfPlan_ = fftw_plan_dft_r2c_1d(n0, &rField[0], &kField[0], flags);
      criPlan_ = fftw_plan_dft_c2r_1d(n0, &kField[0], &rField[0], flags);
      int sign = FFTW_FORWARD;
      ccfPlan_ = fftw_plan_dft_1d(n0, &cFieldIn[0], &cFieldOut[0], 
                                  sign, flags);
      sign = FFTW_BACKWARD;
      cciPlan_ = fftw_plan_dft_1d(n0, &cFieldOut[0], &cFieldIn[0], 
                                  sign, flags);
   }

   template <>
   void FFT<2>::makePlans(RField<2>& rField, RFieldDft<2>& kField,
                          CField<2>& cFieldIn, CField<2>& cFieldOut)
   {
      unsigned int flags = FFTW_ESTIMATE;
      int n0 = meshDimensions_[0];
      int n1 = meshDimensions_[1];
      rcfPlan_ = fftw_plan_dft_r2c_2d(n0, n1, &rField[0], &kField[0], flags);
      criPlan_ = fftw_plan_dft_c2r_2d(n0, n1, &kField[0], &rField[0], flags);
      int sign = FFTW_FORWARD;
      ccfPlan_ = fftw_plan_dft_2d(n0, n1, &cFieldIn[0], &cFieldOut[0], 
                                  sign, flags);
      sign = FFTW_BACKWARD;
      cciPlan_ = fftw_plan_dft_2d(n0, n1, &cFieldOut[0], &cFieldIn[0], 
                                  sign, flags);
   }

   template <>
   void FFT<3>::makePlans(RField<3>& rField, RFieldDft<3>& kField,
                          CField<3>& cFieldIn, CField<3>& cFieldOut)
   {
      unsigned int flags = FFTW_ESTIMATE;
      int n0 = meshDimensions_[0];
      int n1 = meshDimensions_[1];
      int n2 = meshDimensions_[2];
      rcfPlan_ = fftw_plan_dft_r2c_3d(n0, n1, n2, 
                                      &rField[0], &kField[0], flags);
      criPlan_ = fftw_plan_dft_c2r_3d(n0, n1, n2,
                                      &kField[0], &rField[0], flags);
      int sign = FFTW_FORWARD;
      ccfPlan_ = fftw_plan_dft_3d(n0, n1, n2,
                                  &cFieldIn[0], &cFieldOut[0], sign, flags);
      sign = FFTW_BACKWARD;
      cciPlan_ = fftw_plan_dft_3d(n0, n1, n2,
                                  &cFieldIn[0], &cFieldOut[0], sign, flags);
   }

}
}
}
