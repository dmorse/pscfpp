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
      
   double AmStrategyCPU::findResNorm(FieldCPU const & resHist) const 
   {
      const int n = resHist.capacity();
      double normResSq = 0.0;

      for (int i = 0; i < n; i++) {
         normResSq += resHist[i] * resHist[i];
      }

      return sqrt(normResSq);
   }

   double AmStrategyCPU::findResMax(FieldCPU const & resHist) const
   {
      const int n = resHist.capacity();
      double maxRes = 0.0;

      for (int i = 0; i < n; i++) {
         if (fabs(resHist[i]) > maxRes) 
            maxRes = fabs(resHist[i]);
      }

      return maxRes;
   }

   double AmStrategyCPU::computeUDotProd(RingBuffer<FieldCPU> const & resHists, int m) const
   {
      const int n = resHists[0].capacity();
      
      double dotprod = 0.0;
      for(int i = 0; i < n; i++) {
         dotprod += (resHists[m][i] - resHists[m+1][i]) *
                     (resHists[0][i] - resHists[1][i]);
      }

      return dotprod;
   }

   double AmStrategyCPU::computeVDotProd(RingBuffer<FieldCPU> const & resHists, int m) const
   {
      const int n = resHists[0].capacity();
      
      double dotprod = 0.0;
      for(int i = 0; i < n; i++) {
         dotprod += resHists[0][i] * (resHists[m][i] - resHists[m+1][i]);
      }

      return dotprod;
   }

   void AmStrategyCPU::setEqual(FieldCPU& a, FieldCPU const & b) const
   {
      // This seems silly here, but in other implementations it may not be! 
      a = b;
   }

   void AmStrategyCPU::addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & hists, DArray<double> coeffs, int nHist_) const
   {
      int n = trial.capacity();
      for (int i = 0; i < nHist_; i++) {
         for (int j = 0; j < n; j++) {
            trial[j] += coeffs[i] * ( hists[i+1][j] - hists[i][j] );
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