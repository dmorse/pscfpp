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
      
      double findNorm(FieldCPU const & hist) 
      const;

      double findMaxAbs(FieldCPU const & hist) 
      const;

      void updateBasis(RingBuffer<FieldCPU> & basis, RingBuffer<FieldCPU> const & hists)
      const;

      double computeUDotProd(RingBuffer<FieldCPU> const & resBasis, int m, int n) 
      const;

      double computeVDotProd(FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int m) 
      const;

      void updateU(DMatrix<double> & U, RingBuffer<FieldCPU> const & resBasis, int nHist)
      const;

      void updateV(DArray<double> & v, FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int nHist)
      const;

      void setEqual(FieldCPU& a, FieldCPU const & b)
      const;

      void addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & basis, DArray<double> coeffs, int nHist_) 
      const;

      void addPredictedError(FieldCPU& fieldTrial, FieldCPU const & resTrial, double lambda)
      const;

   };

}
}
#endif