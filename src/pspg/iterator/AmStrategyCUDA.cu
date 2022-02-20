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
    : temp_(0)
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
      int nPow2 = n;
      int nExcess = 0;
      if ( (n & (n - 1)) != 0 ) {
         // if not, handle only up to the nearest lower power of two with parallel reduction
         // handle the rest with the CPU
         nPow2 = pow(2,floor(log2(n)));
         nExcess = n-nPow2;
      }

      // Check to verify that private members are allocated
      // Allocate up to n, because other kernels may potentially use these
      // workspace arrays 
      if (!temp_) allocatePrivateMembers(n);

      // Do first reduction step
      int nBlocks = NUMBER_OF_BLOCKS, nThreads=THREADS_PER_BLOCK;
      reductionMaxAbs<<<nBlocks, nThreads, 
                        nThreads*sizeof(cudaReal)>>>(d_temp_, hist.cDField(), nPow2);
      nPow2 = nBlocks;

      // While loop to do further reduction steps
      int itr = 0;
      while (nPow2 > 1) {
         // determine next compute size
         if (nBlocks < nThreads) {
            nThreads = nPow2;
            nBlocks = 1;
         } else {
            nBlocks = nPow2 / nThreads;
         }
         // perform reduction
         reductionMaxAbs<<<nBlocks, nThreads, 
                        nThreads*sizeof(cudaReal)>>>(d_temp_, d_temp_, nPow2);
         nPow2 = nBlocks;

         // track number of iterations
         itr+=1;
         if (itr > 100) {
            UTIL_THROW("Runaway parallel reduction while-loop.");
         }
      }

      cudaReal max;
      cudaMemcpy(&max, d_temp_, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);
      
      // If any left over above the power of two, handle with CPU.
      if (nExcess > 0 ) {
         // copy excess into host memory
         cudaReal* excess = new cudaReal[nExcess];
         cudaMemcpy(excess, hist.cDField() + nPow2, nExcess*sizeof(cudaReal), cudaMemcpyDeviceToHost);
         
         for (int i = 0; i < nExcess; i++) {
            if (fabs(excess[i]) > max) {
               max = fabs(excess[i]);
            }
         }
         delete[] excess;
      }

      return (double)max;

   }

   void AmStrategyCUDA::updateBasis(RingBuffer<FieldCUDA> & basis, RingBuffer<FieldCUDA> const & hists) const
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      FieldCUDA newbasis;
      newbasis.allocate(n);

      pointWiseBinarySubtract<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>
            (hists[0].cDField(),hists[1].cDField(),newbasis.cDField(),n);

      basis.append(newbasis);
   }

   double AmStrategyCUDA::computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int m) const
   {      
      return (double)innerProduct(resBasis[0],resBasis[m]);
   }

   double AmStrategyCUDA::computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m) const
   {
      return (double)innerProduct(resCurrent, resBasis[m]);
   }

   void AmStrategyCUDA::updateU(DMatrix<double> & U, RingBuffer<FieldCUDA> const & resBasis, int nHist) const
   {
      // Update matrix U by shifting elements diagonally
      int maxHist = U.capacity1();
      for (int m = maxHist-1; m > 0; --m) {
         for (int n = maxHist-1; n > 0; --n) {
            U(m,n) = U(m-1,n-1); 
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nHist; ++m) {
         double dotprod = computeUDotProd(resBasis,m);
         U(m,0) = dotprod;
         U(0,m) = dotprod;
      }
   }

   void AmStrategyCUDA::updateV(DArray<double> & v, FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int nHist) const
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
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
               (trial.cDField(), basis[i].cDField(), -1*coeffs[i], trial.capacity());
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
      int n = a.capacity();

      // Check to verify that private members are allocated
      if (!temp_) allocatePrivateMembers(n);

      // Elementwise multiplication of a with b for all elements
      pointWiseBinaryMultiply<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(a.cDField(),b.cDField(),d_temp_,n);

      // Check to see if power of two
      int nPow2 = n;
      int nExcess = 0;
      if ( (n & (n - 1)) != 0 ) {
         // if not, handle only up to the nearest lower power of two with parallel reduction
         // handle the rest with the CPU
         nPow2 = pow(2,floor(log2(n)));
         nExcess = n-nPow2;
      }

      // Parallel reduction (summation over all elements, up to the highest power of two)
      reductionSum<<<NUMBER_OF_BLOCKS/2,THREADS_PER_BLOCK,THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_temp_,d_temp_,nPow2);
      
      // Copy results and sum over output
      cudaMemcpy(temp_, d_temp_, NUMBER_OF_BLOCKS/2 * sizeof(cudaReal), cudaMemcpyDeviceToHost);
      cudaReal sum = 0;
      cudaReal c = 0;
      //use kahan summation to reduce error
      for (int i = 0; i < NUMBER_OF_BLOCKS; ++i) {
         cudaReal y = temp_[i] - c;
         cudaReal t = sum + y;
         c = (t - sum) - y;
         sum = t;  
      }

      // Include data points beyond the lowest power of two
      if (nExcess > 0) {
         // copy excess into host memory
         cudaReal* excess = new cudaReal[nExcess];
         cudaMemcpy(excess, d_temp_ + nPow2, nExcess*sizeof(cudaReal), cudaMemcpyDeviceToHost);
         
         for (int i = 0; i < nExcess; i++) {
            sum += excess[i];
         }
         delete[] excess;
      }
      
      return sum;
   }

}
}



#endif