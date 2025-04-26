/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   // Forward real-to-complex transform, explicit specializations.
   template<>
   void FFT<1>::makePlans()
   {
      int n0 = meshDimensions_[0];
      #ifdef SINGLE_PRECISION
      cufftPlan1d(&rcfPlan_, n0, CUFFT_R2C, 1);
      cufftPlan1d(&criPlan_, n0, CUFFT_C2R, 1);
      cufftPlan1d(&ccPlan_, n0, CUFFT_C2C, 1);
      #else
      cufftPlan1d(&rcfPlan_, n0, CUFFT_D2Z, 1);
      cufftPlan1d(&criPlan_, n0, CUFFT_Z2D, 1);
      cufftPlan1d(&ccPlan_, n0, CUFFT_Z2Z, 1);
      #endif
   }

   template <>
   void FFT<2>::makePlans()
   {
      int n0 = meshDimensions_[0];
      int n1 = meshDimensions_[1];
      #ifdef SINGLE_PRECISION
      cufftPlan2d(&rcfPlan_, n0, n1, CUFFT_R2C);
      cufftPlan2d(&criPlan_, n0, n1, CUFFT_C2R);
      cufftPlan2d(&ccPlan_, n0, n1, CUFFT_C2C);
      #else
      cufftPlan2d(&rcfPlan_, n0, n1, CUFFT_D2Z);
      cufftPlan2d(&criPlan_, n0, n1, CUFFT_Z2D);
      cufftPlan2d(&ccPlan_, n0, n1, CUFFT_Z2Z);
      #endif
   }

   template <>
   void FFT<3>::makePlans()
   {
      int n0 = meshDimensions_[0];
      int n1 = meshDimensions_[1];
      int n2 = meshDimensions_[2];
      #ifdef SINGLE_PRECISION
      cufftPlan3d(&rcfPlan_, n0, n1, n2, CUFFT_R2C);
      cufftPlan3d(&criPlan_, n0, n1, n2, CUFFT_C2R);
      cufftPlan3d(&ccPlan_, n0, n1, n2, CUFFT_C2C);
      #else
      cufftPlan3d(&rcfPlan_, n0, n1, n2, CUFFT_D2Z);
      cufftPlan3d(&criPlan_, n0, n1, n2, CUFFT_Z2D);
      cufftPlan3d(&ccPlan_, n0, n1, n2, CUFFT_Z2Z);
      #endif
   }

   // Explicit instantiation of relevant class instances
   template class FFT<1>;
   template class FFT<2>;
   template class FFT<3>;

}
}
}
