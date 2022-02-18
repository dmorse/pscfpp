#ifndef PSPG_AM_STRATEGY_CUDA_CU
#define PSPG_AM_STRATEGY_CUDA_CU

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmStrategyCUDA.h"
#include <cmath>

namespace Pscf {
namespace Pspg {

   AmStrategyCUDA::AmStrategyCUDA()
   {}

   AmStrategyCUDA::~AmStrategyCUDA()
   {}
      
   double AmStrategyCUDA::findNorm(FieldCUDA const & hist) const 
   {
      const int n = hist.capacity();
      double normResSq = (double)innerProduct(hist, hist);

      return sqrt(normResSq);
   }

   double AmStrategyCUDA::findMaxAbs(FieldCUDA const & hist) const
   {
      // use parallel reduction to find maximum 
      const int n = hist.capacity();
      int outsize = NUMBER_OF_BLOCKS;

      // some loop that changes size of the output as you go, so that parallel reduction is done
      // first for each block, then across the outputs of all the blocks, then across all those outputs,
      // etc!! 
      // use reductionMax 
      // need to modify somehow to find max magnitude, ignoring sign! 
   }

   double AmStrategyCUDA::computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int m) const
   {      
      return (double)innerProduct(resBasis[0],resBasis[m]);
   }

   double AmStrategyCUDA::computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m) const
   {
      return (double)innerProduct(resCurrent, resBasis[m]);
   }

   void AmStrategyCUDA::setEqual(FieldCUDA& a, FieldCUDA const & b) const
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>(a.cDField(), b.cDField(), a.capacity())
   }

   void AmStrategyCUDA::addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & basis, DArray<double> coeffs, int nHist) const
   {
      for (int i = 0; i < nHist; i++) {
         pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> 
               (trial.cDField(), basis[i].cDField(), coeffs[i], trial.capacity());
      }
   }

   void AmStrategyCUDA::addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda) const
   {
      pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> 
         (fieldTrial.cDField(), resTrial.cDField(), lambda, fieldTrial.capacity());
   }

   cudaReal AmStrategyCUDA::innerProduct(FieldCUDA const & a, FieldCUDA const & b) const
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      int size = a.capacity();

      switch(THREADS_PER_BLOCK) {
      case 512:
         deviceInnerProduct<512><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 256:
         deviceInnerProduct<256><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 128:
         deviceInnerProduct<128><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 64:
         deviceInnerProduct<64><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 32:
         deviceInnerProduct<32><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 16:
         deviceInnerProduct<16><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 8:
         deviceInnerProduct<8><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 4:
         deviceInnerProduct<4><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 2:
         deviceInnerProduct<2><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      case 1:
         deviceInnerProduct<1><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
         break;
      }
      cudaMemcpy(temp_, d_temp_, NUMBER_OF_BLOCKS * sizeof(cudaReal), cudaMemcpyDeviceToHost);
      cudaReal final = 0;
      cudaReal c = 0;
      //use kahan summation to reduce error
      for (int i = 0; i < NUMBER_OF_BLOCKS; ++i) {
         cudaReal y = temp_[i] - c;
         cudaReal t = final + y;
         c = (t - final) - y;
         final = t;  
      }
      return final;
   }

}
}



#endif