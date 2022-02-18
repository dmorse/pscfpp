#ifndef PSPG_AM_STRATEGY_CUDA_H
#define PSPG_AM_STRATEGY_CUDA_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmStrategy.h>
#include <util/global.h>
#include <pspg/field/DField.h>
#include <pspg/math/GpuResources.h>
#include <util/containers/RingBuffer.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   typedef DField<cudaReal> FieldCUDA;

   class AmStrategyCUDA : public AmStrategy<FieldCUDA>
   {
   public:
         /// Constructor
      AmStrategyCUDA();

      /// Destructor
      ~AmStrategyCUDA(); 
      
      double findResNorm(FieldCUDA const & resHist) 
      const;

      double findResMax(FieldCUDA const & resHist) 
      const;

      double computeUDotProd(RingBuffer<FieldCUDA> const & resHists, int m) 
      const;

      double computeVDotProd(RingBuffer<FieldCUDA> const & resHists, int m) 
      const;

      void setEqual(FieldCUDA& a, FieldCUDA const & b)
      const;

      void addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & hists, DArray<double> coeffs, int nHist_) 
      const;

      void addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda)
      const;

   };

}
}
#endif