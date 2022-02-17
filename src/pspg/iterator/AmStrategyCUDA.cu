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
      
   double AmStrategyCUDA::findResNorm(FieldCUDA const & resHist) const 
   {

   }

   double AmStrategyCUDA::findResMax(FieldCUDA const & resHist) const
   {

   }

   double AmStrategyCUDA::computeUDotProd(RingBuffer<FieldCUDA> const & resHists, int m) const
   {

   }

   double AmStrategyCUDA::computeVDotProd(RingBuffer<FieldCUDA> const & resHists, int m) const
   {

   }

   void AmStrategyCUDA::setEqual(FieldCUDA& a, FieldCUDA const & b) const
   {

   }

   void AmStrategyCUDA::addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & hists, DArray<double> coeffs, int nHist_) const
   {

   }

   void AmStrategyCUDA::addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda) const
   {

   }

}
}
#endif