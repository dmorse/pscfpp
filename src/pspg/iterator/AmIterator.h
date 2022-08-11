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
      *   
      * \param system parent system object
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

      using AmIteratorTmpl<Iterator<D>,FieldCUDA>::setup;
      using AmIteratorTmpl<Iterator<D>,FieldCUDA>::solve;
      using Iterator<D>::isFlexible;

   protected:

      using ParamComposite::readOptional;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;

   private:

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

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
      * \param basis RingBuffer of basis vectors.
      * \param hists RingBuffer of previous vectors.
      */
      void updateBasis(RingBuffer<FieldCUDA> & basis, 
                       RingBuffer<FieldCUDA> const & hists);

      /**
      * Compute the dot product for an element of the U matrix.
      * 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, 
                             int m, int n);

      /**
      * Compute the dot product for an element of the v vector.
      * 
      * \param resCurrent current residual vector
      * \param resBasis RingBuffer of residual basis vectors
      * \param m row index for element of the v vector
      */
      double computeVDotProd(FieldCUDA const & resCurrent, 
                             RingBuffer<FieldCUDA> const & resBasis, 
                             int m);

      /**
      * Update the U matrix.
      * 
      * \param U U matrix (dot products of residual basis vectors)
      * \param resBasis RingBuffer residual basis vectors
      * \param nHist current number of previous states
      */
      void updateU(DMatrix<double> & U, 
                   RingBuffer<FieldCUDA> const & resBasis, 
                   int nHist);

      /**
      * Update the v vector.
      * 
      * \param v v vector
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void updateV(DArray<double> & v, 
                   FieldCUDA const & resCurrent, 
                   RingBuffer<FieldCUDA> const & resBasis, 
                   int nHist);

      /**
      * Set a vector equal to another. 
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(FieldCUDA& a, FieldCUDA const & b);

      /**
      * Compute trial field so as to minimize L2 norm of residual.
      * 
      * \param trial resulting trial field (output)
      * \param basis RingBuffer of residual basis vectors.
      * \param coeffs coefficients of basis vectors
      * \param nHist number of prior states stored
      */
      void addHistories(FieldCUDA& trial, 
                        RingBuffer<FieldCUDA> const & basis, 
                        DArray<double> coeffs, int nHist);

      /**
      * Add predicted error to the trial field.
      * 
      * \param fieldTrial trial field (input/output)
      * \param resTrial predicted error for current trial field
      * \param lambda Anderson-Mixing mixing parameter 
      */
      void addPredictedError(FieldCUDA& fieldTrial, 
                             FieldCUDA const & resTrial, 
                             double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
     
      /** 
      * Compute the number of elements in field or residual.
      */
      int nElements();

      /*
      * Get the current state of the system.
      *
      * \param curr current field vector (output)
      */
      void getCurrent(FieldCUDA& curr);

      /**
      * Solve MDE for current state of system.
      */
      void evaluate();

      /**
      * Gets the residual vector from system.
      *  
      * \param curr current residual vector (output)
      */
      void getResidual(FieldCUDA& resid);

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess trial field configuration
      */
      void update(FieldCUDA& newGuess);

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog();

      // --- Private member functions specific to this implementation --- 
      
      cudaReal findAverage(cudaReal * const field, int n);

   };

} // namespace Pspg
} // namespace Pscf
#endif
