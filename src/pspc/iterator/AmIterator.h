#ifndef PSPC_AM_ITERATOR_H
#define PSPC_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspc implementation of the Anderson Mixing iterator.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class AmIterator : public AmIteratorTmpl<Iterator<D>, FieldCPU>
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this iterator.
      */
      AmIterator(System<D>& system);

      /**
      * Destructor.
      */
      ~AmIterator();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      using AmIteratorTmpl<Iterator<D>,FieldCPU>::setup;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::solve;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::maskField;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::externalFields;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::externalField;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::hasMask;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::hasExternalField;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::setClassName;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;

   protected:
   
      using ParamComposite::readOptional;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::flexibleParams_;

   private:

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;
      
      /**
      * Find norm of a residual vector.
      */
      double findNorm(FieldCPU const & hist);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double findMaxAbs(FieldCPU const & hist);

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer object storing the list of residual or field basis vectors.
      * \param hists RingBuffer object storing the histories of residual or field vectors.
      */
      void updateBasis(RingBuffer<FieldCPU> & basis, RingBuffer<FieldCPU> const & hists);

      /**
      * Compute the dot product for an element of the U matrix.
      * 
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<FieldCPU> const & resBasis, int m, int n);

      /**
      * Compute the dot product for an element of the v vector.
      * 
      * \param resCurrent the residual vector calculated at the present iteration step
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param m row of the v vector
      */
      double computeVDotProd(FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int m);

      /**
      * Compute the series of necessary dot products and update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      void updateU(DMatrix<double> & U, RingBuffer<FieldCPU> const & resBasis, int nHist);

      /**
      * Compute the series of necessary dot products and update the v vector.
      * 
      * \param v v vector
      * \param resCurrent the residual vector calculated at the present iteration step
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      void updateV(DArray<double> & v, FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int nHist);

      /**
      * Set a field equal to another. Essentially a = b, but potentially more complex
      * in certain implementations of the AmIterator.
      * 
      * \param a the field to be set
      * \param b the field for it to be set to
      */
      void setEqual(FieldCPU& a, FieldCPU const & b);

      /**
      * Mix histories, scaled by their respective coefficients, into the trial field.
      * 
      * \param trial object for calculation results to be stored in.
      * \param basis list of history basis vectors.
      * \param coeffs list of coefficients for each history.
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & basis, DArray<double> coeffs, int nHist);

      /**
      * Add predicted error into the field trial guess to attempt to correct for it.
      * 
      * \param fieldTrial field for calculation results to be stored in.
      * \param resTrial predicted error for current mixing of histories.
      * \param lambda Anderson-Mixing parameter for mixing in histories
      */
      void addPredictedError(FieldCPU& fieldTrial, FieldCPU const & resTrial, double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated
      int nElements();

      /// Gets the current state of the system and processes it into 
      /// format needed by AmIteratorTmpl.
      void getCurrent(FieldCPU& curr);

      /// Runs calculation to evaluate function for fixed point.
      void evaluate();

      /// Gets residual values from system
      void getResidual(FieldCPU& resid);

      /// Updates the system with a passed in state of the iterator.
      void update(FieldCPU& newGuess);

      /// Outputs relevant system details to the iteration log
      void outputToLog();

   };

} // namespace Pspc
} // namespace Pscf
#endif
