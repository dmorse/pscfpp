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
      using Iterator<D>::sys_;

   private:

      /**
      * Find norm of a residual vector.
      */
      double findNorm(FieldCUDA const & hist);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double findMaxAbs(FieldCUDA const & hist);

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer object storing the list of residual or field basis vectors.
      * \param hists RingBuffer object storing the histories of residual or field vectors.
      */
      void updateBasis(RingBuffer<FieldCUDA> & basis, RingBuffer<FieldCUDA> const & hists);

      /**
      * Compute the dot product for an element of the U matrix.
      * 
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int m, int n);

      /**
      * Compute the dot product for an element of the v vector.
      * 
      * \param resCurrent the residual vector calculated at the present iteration step
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param m row of the v vector
      */
      double computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m);

      /**
      * Compute the series of necessary dot products and update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      void updateU(DMatrix<double> & U, RingBuffer<FieldCUDA> const & resBasis, int nHist);

      /**
      * Compute the series of necessary dot products and update the v vector.
      * 
      * \param v v vector
      * \param resCurrent the residual vector calculated at the present iteration step
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      void updateV(DArray<double> & v, FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int nHist);

      /**
      * Set a field equal to another. Essentially a = b, but potentially more complex
      * in certain implementations of the AmIterator.
      * 
      * \param a the field to be set
      * \param b the field for it to be set to
      */
      void setEqual(FieldCUDA& a, FieldCUDA const & b);

      /**
      * Mix histories, scaled by their respective coefficients, into the trial field.
      * 
      * \param trial object for calculation results to be stored in.
      * \param basis list of history basis vectors.
      * \param coeffs list of coefficients for each history.
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & basis, DArray<double> coeffs, int nHist);

      /**
      * Add predicted error into the field trial guess to attempt to correct for it.
      * 
      * \param fieldTrial field for calculation results to be stored in.
      * \param resTrial predicted error for current mixing of histories.
      * \param lambda Anderson-Mixing parameter for mixing in histories
      */
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
