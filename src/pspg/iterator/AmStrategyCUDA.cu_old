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
      double normResSq = (double)gpuInnerProduct(hist.cDField(), hist.cDField(), n);

      return sqrt(normResSq);
   }

   double AmStrategyCUDA::findMaxAbs(FieldCUDA const & hist) const
   {
      // use parallel reduction to find maximum.

      // number of data points, each step of the way.
      int n = hist.capacity();
      cudaReal max = gpuMaxAbs(hist.cDField(), n);

      return (double)max;

   }

   void AmStrategyCUDA::updateBasis(RingBuffer<FieldCUDA> & basis, RingBuffer<FieldCUDA> const & hists) const
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      // Clear current basis vectors, determine number to compute
      basis.clear();
      const int nBasisVec = hists.size() - 1;
      const int n = hists[0].capacity();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Holding area for new basis vectors
      FieldCUDA newbasis;
      newbasis.allocate(n);

      // Compute each basis vector
      for (int i = 0; i < nBasisVec; i++)
      {
         pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            (hists[0].cDField(),hists[i+1].cDField(),newbasis.cDField(),n);
            basis.append(newbasis);
      }
      
   }

   double AmStrategyCUDA::computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int n, int m) const
   {      
      return (double)gpuInnerProduct(resBasis[n].cDField(),resBasis[m].cDField(), resBasis[n].capacity());
   }

   double AmStrategyCUDA::computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m) const
   {
      return (double)gpuInnerProduct(resCurrent.cDField(), resBasis[m].cDField(), resCurrent.capacity());
   }

   void AmStrategyCUDA::updateU(DMatrix<double> & U, RingBuffer<FieldCUDA> const & resBasis, int nHist) const
   {
      // Reset U values to 0
      int maxHist = U.capacity1();
      for (int m = 0; m < maxHist; m++) {
         for (int n = 0; n < maxHist; n++) {
            U(m,n) = 0;
         }
      }

      // Compute U matrix's new rows and columns and col 0
      for (int m = 0; m < nHist; ++m) {
         for (int n = m; n < nHist; ++n) {
            double dotprod = computeUDotProd(resBasis,m,n);
            U(m,n) = dotprod;
            U(n,m) = dotprod;
         }
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
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cDField(), b.cDField(), a.capacity());
   }

   void AmStrategyCUDA::addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & basis, DArray<double> coeffs, int nHist) const
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(trial.capacity(), nBlocks, nThreads);

      for (int i = 0; i < nHist; i++) {
         pointWiseAddScale<<<nBlocks, nThreads>>>
               (trial.cDField(), basis[i].cDField(), -1*coeffs[i], trial.capacity());
      }
   }

   void AmStrategyCUDA::addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda) const
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);

      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cDField(), resTrial.cDField(), lambda, fieldTrial.capacity());
   }

}
}



#endif