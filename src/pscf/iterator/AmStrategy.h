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
   public:

      /// Constructor
      AmStrategy();

      /// Destructor
      ~AmStrategy(); 
      
      virtual double findResNorm(T const & resHist) const = 0;
      virtual double findResMax(T const & resHist) const = 0;
      virtual double computeUDotProd(RingBuffer<T> const & resHists, int m) const = 0;
      virtual double computeVDotProd(RingBuffer<T> const & resHists, int m) const = 0;
      virtual void setEqual(T& a, T const & b) const = 0;
      virtual void addHistories(T& trial, RingBuffer<T> const & hists, DArray<double> coeffs, int nHist_) const = 0;
      virtual void addPredictedError(T& fieldTrial, T const & resTrial, double lambda) const = 0;

   };

}
#endif