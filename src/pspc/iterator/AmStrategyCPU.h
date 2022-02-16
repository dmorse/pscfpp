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

   class AmStrategyCPU : public AmStrategy<FieldCPU>
   {
   public:
         /// Constructor
      AmStrategyCPU();

      /// Destructor
      ~AmStrategyCPU(); 
      
      double findResNorm(FieldCPU const & resHist) 
      const;

      double findResMax(FieldCPU const & resHist) 
      const;

      double computeUDotProd(RingBuffer<FieldCPU> const & resHists, int m) 
      const;

      double computeVDotProd(RingBuffer<FieldCPU> const & resHists, int m) 
      const;

      void setEqual(FieldCPU& a, FieldCPU const & b)
      const;

      void addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & hists, DArray<double> coeffs, int nHist_) 
      const;

      void addPredictedError(FieldCPU& fieldTrial, FieldCPU const & resTrial, double lambda)
      const;

   };

}
}
#endif