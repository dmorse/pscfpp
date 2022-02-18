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
   {
      if (temp_) {
         delete[] temp_;
         cudaFree(d_temp_);
      }
   }
      
   double AmStrategyCUDA::findNorm(FieldCUDA const & hist) const 
   {
      const int n = hist.capacity();
      double normResSq = (double)innerProduct(hist, hist);

      return sqrt(normResSq);
   }

   double AmStrategyCUDA::findMaxAbs(FieldCUDA const & hist) const
   {
      // use parallel reduction to find maximum.

      // number of data points, each step of the way.
      int n = hist.capacity();

      // CHECK IF POWER OF TWO. 
      // bitwise trick from http://www.graphics.stanford.edu/~seander/bithacks.html
      if ( (n & (n - 1)) != 0 )
         UTIL_THROW("Dataset size must be a power of two.");

      // Check to verify that private members are allocated
      if (!temp_) allocatePrivateMembers(n);

      // Do first reduction step
      int nBlocks = NUMBER_OF_BLOCKS, nThreads=THREADS_PER_BLOCK;
      reductionMaxAbs<<<nBlocks, nThreads, 
                        nThreads*sizeof(cudaReal)>>>(d_temp_, hist.cDField(), n);
      n = nBlocks;

      // While loop to do further reduction steps
      int itr = 0;
      while (n > 1) {
         // determine next compute size
         if (nBlocks < nThreads) {
            nThreads = n;
            nBlocks = 1;
         } else {
            nBlocks = n / nThreads;
         }
         // perform reduction
         reductionMaxAbs<<<nBlocks, nThreads, 
                        nThreads*sizeof(cudaReal)>>>(d_temp_, d_temp_, n);
         n = nBlocks;

         // track number of iterations
         itr+=1;
         if (itr > 100) {
            UTIL_THROW("Runaway parallel reduction while-loop.");
         }
      }

      cudaReal max;
      cudaMemcpy(&max, d_temp_, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);
      return (double)max;

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
      assignReal<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>(a.cDField(), b.cDField(), a.capacity());
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

   // --- Private member functions that are specific to this implementation --- 

   void AmStrategyCUDA::allocatePrivateMembers(int n) const
   {
      temp_ = new cudaReal[n];
      cudaMalloc((void**) &d_temp_, n*sizeof(cudaReal));
   }

   cudaReal AmStrategyCUDA::innerProduct(FieldCUDA const & a, FieldCUDA const & b) const
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      int size = a.capacity();

      // Check to verify that private members are allocated
      if (!temp_) allocatePrivateMembers(size);

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