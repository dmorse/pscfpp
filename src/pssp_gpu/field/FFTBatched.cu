/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.h"

namespace Pscf {
namespace Pssp_gpu {

   using namespace Util;

   // Forward transform, explicit specializations.
   //use local mesh dimensions later
   /*template<>
   void FFTBatched<1>::makePlans(RDField<1>& rField, RDFieldDft<1>& kField)
   {
      cufftPlan1d(&fPlan_, rField.capacity(), CUFFT_R2C, 1);
      cufftPlan1d(&iPlan_, rField.capacity(), CUFFT_C2R, 1);
   }

   template <>
   void FFTBatched<2>::makePlans(RDField<2>& rField, RDFieldDft<2>& kField)
   {
      cufftPlan2d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_R2C);
      cufftPlan2d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1], CUFFT_C2R);

      }*/

   template <>
   void FFTBatched<3>::makePlans(RDField<3>& rField, RDFieldDft<3>& kField)
   {
      //std::cout<<"rfielddim2 "<<rField.meshDimensions()[2]<<std::endl;
      //cufftPlan3d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
      //      rField.meshDimensions()[2], CUFFT_R2C);
      //cufftPlan3d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
      //     rField.meshDimensions()[2], CUFFT_C2R);
      int n[3];
      n[0] = rField.meshDimensions()[0];
      n[1] = rField.meshDimensions()[1];
      n[2] = rField.meshDimensions()[2];
      if(cufftPlanMany(&fPlan_, 3, n, //plan, rank, n
                       NULL, 1, n[0]*n[1]*n[2], //inembed, istride, idist
                       NULL, 1, n[0]*n[1]*(n[2]/2 + 1), //onembed, ostride, odist
                       CUFFT_R2C, 2) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if(cufftPlanMany(&iPlan_, 3, n, //plan, rank, n
                       NULL, 1, n[0]*n[1]*(n[2]/2 + 1), //inembed, istride, idist
                       NULL, 1, n[0]*n[1]*n[2], //onembed, ostride, odist
                       CUFFT_C2R, 2) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      
   }

   template <>
   void FFTBatched<3>::makePlans(const IntVec<3>& rDim, const IntVec<3>& kDim, int batchSize)
   {
      //std::cout<<"rfielddim2 "<<rField.meshDimensions()[2]<<std::endl;
      //cufftPlan3d(&fPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
      //      rField.meshDimensions()[2], CUFFT_R2C);
      //cufftPlan3d(&iPlan_, rField.meshDimensions()[0], rField.meshDimensions()[1],
      //     rField.meshDimensions()[2], CUFFT_C2R);
      int n[3];
      n[0] = rDim[0];
      n[1] = rDim[1];
      n[2] = rDim[2];
      if(cufftPlanMany(&fPlan_, 3, n, //plan, rank, n
                       NULL, 1, n[0]*n[1]*n[2], //inembed, istride, idist
                       NULL, 1, n[0]*n[1]*(n[2]/2 + 1), //onembed, ostride, odist
                       CUFFT_R2C, batchSize) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if(cufftPlanMany(&iPlan_, 3, n, //plan, rank, n
                       NULL, 1, n[0]*n[1]*(n[2]/2 + 1), //inembed, istride, idist
                       NULL, 1, n[0]*n[1]*n[2], //onembed, ostride, odist
                       CUFFT_C2R, batchSize) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      
   }
}
}
