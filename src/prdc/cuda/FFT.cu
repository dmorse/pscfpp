/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   // Forward real-to-complex transform, explicit specializations.
   template<>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField)
   {
      #ifdef SINGLE_PRECISION
      cufftPlan1d(&rcfPlan_, rField.capacity(), CUFFT_R2C, 1);
      cufftPlan1d(&criPlan_, rField.capacity(), CUFFT_C2R, 1);
      #else
      cufftPlan1d(&rcfPlan_, rField.capacity(), CUFFT_D2Z, 1);
      cufftPlan1d(&criPlan_, rField.capacity(), CUFFT_Z2D, 1);
      #endif
   }

   template <>
   void FFT<2>::makePlans(RField<2>& rField, RFieldDft<2>& kField)
   {
      #ifdef SINGLE_PRECISION
      cufftPlan2d(&rcfPlan_, 
                  rField.meshDimensions()[0], rField.meshDimensions()[1], 
                  CUFFT_R2C);
      cufftPlan2d(&criPlan_, 
                  rField.meshDimensions()[0], rField.meshDimensions()[1], 
                  CUFFT_C2R);
      #else
      cufftPlan2d(&rcfPlan_, 
                  rField.meshDimensions()[0], rField.meshDimensions()[1], 
                  CUFFT_D2Z);
      cufftPlan2d(&criPlan_, 
                  rField.meshDimensions()[0], rField.meshDimensions()[1], 
                  CUFFT_Z2D);
      #endif
   }

   template <>
   void FFT<3>::makePlans(RField<3>& rField, RFieldDft<3>& kField)
   {
      #ifdef SINGLE_PRECISION
      cufftPlan3d(&rcfPlan_, 
                  rField.meshDimensions()[0], 
                  rField.meshDimensions()[1],
                  rField.meshDimensions()[2], CUFFT_R2C);
      cufftPlan3d(&criPlan_, 
                  rField.meshDimensions()[0], 
                  rField.meshDimensions()[1],
                  rField.meshDimensions()[2], CUFFT_C2R);
      #else
      cufftPlan3d(&rcfPlan_, 
                  rField.meshDimensions()[0], 
                  rField.meshDimensions()[1],
                  rField.meshDimensions()[2], CUFFT_D2Z);
      cufftPlan3d(&criPlan_, 
                  rField.meshDimensions()[0], 
                  rField.meshDimensions()[1],
                  rField.meshDimensions()[2], CUFFT_Z2D);
      #endif
   }

   // Explicit instantiation of relevant class instances
   template class FFT<1>;
   template class FFT<2>;
   template class FFT<3>;

}
}
}
