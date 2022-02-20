#ifndef PSPC_AM_STRATEGY_CPU_CPP
#define PSPC_AM_STRATEGY_CPU_CPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmStrategyCPU.h"
#include <cmath>

namespace Pscf {
namespace Pspc {

   AmStrategyCPU::AmStrategyCPU()
   {}

   AmStrategyCPU::~AmStrategyCPU()
   {}
      
   double AmStrategyCPU::findNorm(FieldCPU const & hist) const 
   {
      const int n = hist.capacity();
      double normResSq = 0.0;

      for (int i = 0; i < n; i++) {
         normResSq += hist[i] * hist[i];
      }

      return sqrt(normResSq);
   }

   double AmStrategyCPU::findMaxAbs(FieldCPU const & hist) const
   {
      const int n = hist.capacity();
      double maxRes = 0.0;

      for (int i = 0; i < n; i++) {
         if (fabs(hist[i]) > maxRes) 
            maxRes = fabs(hist[i]);
      }

      return maxRes;
   }

   void AmStrategyCPU::updateBasis(RingBuffer<FieldCPU> & basis, RingBuffer<FieldCPU> const & hists) const
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      FieldCPU newbasis;
      newbasis.allocate(n);

      for (int i = 0; i < n; i++) {
         newbasis[i] = hists[0][i] - hists[1][i]; // sequential histories basis vectors
      }

      basis.append(newbasis);
   }

   double AmStrategyCPU::computeUDotProd(RingBuffer<FieldCPU> const & resBasis, int m, int n) const
   {
      const int n = resBasis[0].capacity();
      
      double dotprod = 0.0;
      for(int i = 0; i < n; i++) {
         dotprod += resBasis[m][i] * resBasis[n][i];
      }

      return dotprod;
   }

   double AmStrategyCPU::computeVDotProd(FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int m) const
   {
      const int n = resBasis[0].capacity();
      
      double dotprod = 0.0;
      for(int i = 0; i < n; i++) {
         dotprod += resCurrent[i] * resBasis[m][i];
      }

      return dotprod;
   }

   void AmStrategyCPU::updateU(DMatrix<double> & U, RingBuffer<FieldCPU> const & resBasis, int nHist) const
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
         double dotprod = computeUDotProd(resBasis,0,m);
         U(m,0) = dotprod;
         U(0,m) = dotprod;
      }
   }

   void AmStrategyCPU::updateV(DArray<double> & v, FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int nHist) const
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }

   void AmStrategyCPU::setEqual(FieldCPU& a, FieldCPU const & b) const
   {
      // This seems silly here, but in other implementations it may not be! Goal: no explicit math in AmIterator.
      a = b;
   }

   void AmStrategyCPU::addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & basis, DArray<double> coeffs, int nHist) const
   {
      int n = trial.capacity();
      for (int i = 0; i < nHist; i++) {
         for (int j = 0; j < n; j++) {
            // Not clear on the origin of the -1 factor
            trial[j] += coeffs[i] * -1 * basis[i][j]; 
         }
      }
   }

   void AmStrategyCPU::addPredictedError(FieldCPU& fieldTrial, FieldCPU const & resTrial, double lambda) const
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
   }

}
}
#endif