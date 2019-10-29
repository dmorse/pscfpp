/*
* PSCF++ Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFT.h"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Forward transform, explicit specializations.
   //use local mesh dimensions later
   template<>
   void FFT<1>::makePlans(RDField<1>& rField, RDFieldDft<1>& kField)
   {
      cufftPlan1d(&fPlan_, rField.capacity(), CUFFT_R2C, 1);
      cufftPlan1d(&iPlan_, rField.capacity(), CUFFT_C2R, 1);
   }

   template <>
   void FFT<2>::makePlans(RDField<2>& rField, RDFieldDft<2>& kField)
   {
      cufftPlan2d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_R2C);
      cufftPlan2d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_C2R);

   }

   template <>
   void FFT<3>::makePlans(RDField<3>& rField, RDFieldDft<3>& kField)
   {
      //std::cout<<"rfielddim2 "<<rField.meshDimensions()[2]<<std::endl;
      cufftPlan3d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
            rField.meshDimensions()[2], CUFFT_R2C);
      cufftPlan3d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
            rField.meshDimensions()[2], CUFFT_C2R);
   }
}
}
