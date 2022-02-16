#ifndef PSCF_AM_STRATEGY_H
#define PSCF_AM_STRATEGY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/containers/DArray.h>
#include <util/containers/RingBuffer.h>

namespace Pscf {

   using namespace Util;

   template <typename T>
   class AmStrategy 
   {
      /// Constructor
      AmStrategy();

      /// Destructor
      ~AmStrategy(); 

   public:
      
      virtual double findResNorm(T& resHist) = 0;
      virtual double findResMax(T& resHist) = 0;
      virtual double computeUDotProd(RingBuffer<T>& resHists, int m) = 0;
      virtual double computeVDotProd(RingBuffer<T>& resHists, int m) = 0;
      virtual void setEqual(T& a, T& b) = 0;
      virtual void addHistories(T& trial, RingBuffer<T>& hists, DArray<double> coeffs, int nHist_) = 0;
      virtual void addPredictedError(T& trial, T& resTrial, double lambda) = 0;

   };

}
#endif