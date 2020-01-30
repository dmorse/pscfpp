/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.tpp"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Forward transform, explicit specializations.
   //use local mesh dimensions later
   template<>
   void FFT<1>::makePlans(RDField<1>& rField, RDFieldDft<1>& kField)
   {
#ifdef SINGLE_PRECISION
      cufftPlan1d(&fPlan_, rField.capacity(), CUFFT_R2C, 1);
      cufftPlan1d(&iPlan_, rField.capacity(), CUFFT_C2R, 1);
#else
      cufftPlan1d(&fPlan_, rField.capacity(), CUFFT_D2Z, 1);
      cufftPlan1d(&iPlan_, rField.capacity(), CUFFT_Z2D, 1);
#endif
   }

   template <>
   void FFT<2>::makePlans(RDField<2>& rField, RDFieldDft<2>& kField)
   {
#ifdef SINGLE_PRECISION
      cufftPlan2d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_R2C);
      cufftPlan2d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_C2R);
#else
      cufftPlan2d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_D2Z);
      cufftPlan2d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_Z2D);
#endif

   }

   template <>
   void FFT<3>::makePlans(RDField<3>& rField, RDFieldDft<3>& kField)
   {
#ifdef SINGLE_PRECISION
      cufftPlan3d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
            rField.meshDimensions()[2], CUFFT_R2C);
      cufftPlan3d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
            rField.meshDimensions()[2], CUFFT_C2R);
#else
      cufftPlan3d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
            rField.meshDimensions()[2], CUFFT_D2Z);
      cufftPlan3d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
            rField.meshDimensions()[2], CUFFT_Z2D);
#endif
   }

   // Explicit instantiation of relevant class instances
   template class FFT<1>;
   template class FFT<2>;
   template class FFT<3>;

}
}
