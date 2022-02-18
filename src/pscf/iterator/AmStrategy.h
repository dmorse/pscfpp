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
      AmStrategy() {};

      /// Destructor
      virtual ~AmStrategy() {}; 
      
      /// Find the norm of the residual vector.
      virtual double findNorm(T const & hist) const = 0;

      /// Find the element of the residual vector with the maximum magnitude.
      virtual double findMaxAbs(T const & hist) const = 0;

      /// Update the list of residual basis vectors used for combining histories.
      virtual void updateBasis(RingBuffer<T> & basis, RingBuffer<T> const & hists) const = 0;
      
      /// Compute the dot product for constructing the U matrix. 
      virtual double computeUDotProd(RingBuffer<T> const & resBasis, int m) const = 0;
      
      /// Compute the dot product for constructing the v vector. 
      virtual double computeVDotProd(T const & resCurrent, RingBuffer<T> const & resBasis, int m) const = 0;
      
      /// Set two things equal to each other.
      virtual void setEqual(T& a, T const & b) const = 0;

      /// Mix histories, scaled by their respective coefficients, into the trial field.
      virtual void addHistories(T& trial, RingBuffer<T> const & basis, DArray<double> coeffs, int nHist) const = 0;

      /// Add predicted error into the trial guess to approximately correct for it.
      virtual void addPredictedError(T& fieldTrial, T const & resTrial, double lambda) const = 0;

   };

}
#endif