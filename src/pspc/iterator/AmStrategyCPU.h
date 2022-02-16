#ifndef PSPC_AM_STRATEGY_CPU_H
#define PSPC_AM_STRATEGY_CPU_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmStrategy.h>
#include <util/global.h>
#include <util/containers/DArray.h>
#include <util/containers/RingBuffer.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   typedef DArray<double> FieldCPU;

   class AmStrategyCPU : AmStrategy<FieldCPU>
   {
      /// Constructor
      AmStrategyCPU();

      /// Destructor
      ~AmStrategyGPU(); 

   public:
      
      virtual double findResNorm(FieldCPU& resHist);

      virtual double findResMax(FieldCPU& resHist);

      virtual double computeUDotProd(RingBuffer<FieldCPU>& resHists, int m);

      virtual double computeVDotProd(RingBuffer<FieldCPU>& resHists, int m);

      virtual void setEqual(FieldCPU& a, FieldCPU& b);

      virtual void addHistories(FieldCPU& trial, RingBuffer<FieldCPU>& hists, DArray<double> coeffs, int nHist_);

      virtual void addPredictedError(FieldCPU& trial, FieldCPU& resTrial, double lambda);


   };

}
}
#endif