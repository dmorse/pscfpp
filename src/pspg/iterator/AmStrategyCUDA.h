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
      
      double findNorm(FieldCUDA const & hist) 
      const;

      double findRelNorm(FieldCUDA const & hist) 
      const;

      double findMaxAbs(FieldCUDA const & hist) 
      const;

      void updateBasis(RingBuffer<FieldCUDA> & basis, RingBuffer<FieldCUDA> const & hists)
      const;

      double computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int m, int n) 
      const;

      double computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m) 
      const;

      void updateU(DMatrix<double> & U, RingBuffer<FieldCUDA> const & resBasis, int nHist)
      const;

      void updateV(DArray<double> & v, FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int nHist)
      const;

      void setEqual(FieldCUDA& a, FieldCUDA const & b)
      const;

      void addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & basis, DArray<double> coeffs, int nHist) 
      const;

      void addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda)
      const;
   
   private:

      void allocatePrivateMembers(int n) const;

      // workspace members
      mutable cudaReal* d_temp_;

      mutable cudaReal* temp_;

   };

}
}
#endif