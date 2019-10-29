/*
* PSCF++ Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.h"

namespace Pscf {
namespace Pspg {

   using namespace Util;

#if 0
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
#endif

#if 0
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
#endif

   template class FFTBatched<1>;
   template class FFTBatched<2>;
   template class FFTBatched<3>;

}
}
