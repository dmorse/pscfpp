#ifndef PSPG_AM_ITERATOR_H
#define PSPG_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspg implementation of the Anderson Mixing iterator.
   *
   * \ingroup Pspg_Iterator_Module
   */
   template <int D>
   class AmIterator : public AmIteratorTmpl<Iterator<D>, FieldCUDA>
   {

   public:
      /**
      * Constructor.
      */
      AmIterator(System<D>& system);

      ~AmIterator();

      using AmIteratorTmpl<Iterator<D>,FieldCUDA>::setup;
      using AmIteratorTmpl<Iterator<D>,FieldCUDA>::solve;
      using AmIteratorTmpl<Iterator<D>,FieldCUDA>::readParameters;

   private:

      System<D>* sys_;

      double findNorm(FieldCUDA const & hist);

      double findRelNorm(FieldCUDA const & hist);

      double findMaxAbs(FieldCUDA const & hist);

      void updateBasis(RingBuffer<FieldCUDA> & basis, RingBuffer<FieldCUDA> const & hists);

      double computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int m, int n);

      double computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m);

      void updateU(DMatrix<double> & U, RingBuffer<FieldCUDA> const & resBasis, int nHist);

      void updateV(DArray<double> & v, FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int nHist);

      void setEqual(FieldCUDA& a, FieldCUDA const & b);

      void addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & basis, DArray<double> coeffs, int nHist);

      void addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated
      int nElements();

      /// Gets a reference to the current state of the system
      void getCurrent(FieldCUDA& curr);

      /// Runs calculation to evaluate function for fixed point.
      void evaluate();

      /// Gets residual values from system
      void getResidual(FieldCUDA& resid);

      /// Updates the system with a passed in state of the iterator.
      void update(FieldCUDA& newGuess);

      /// Outputs relevant system details to the iteration log
      void outputToLog();

      // --- Private member functions that are specific to this implementation --- 
      cudaReal findAverage(cudaReal * const field, int n);

   };

} // namespace Pspg
} // namespace Pscf
#endif
