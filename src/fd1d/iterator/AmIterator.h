#ifndef FD1D_AM_ITERATOR_H
#define FD1D_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <fd1d/solvers/Mixture.h>
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Fd1d
{

   class System;

   using namespace Util;

   /**
   * Anderson-Mixing iterator.
   *
   * \ingroup Fd1d_Iterator_Module
   */
   class AmIterator : public AmIteratorTmpl<Iterator, DArray<double> >
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this iterator.
      */
      AmIterator(System& system);

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
      
      // Inherited public member functions
      using AmIteratorTmpl<Iterator,DArray<double> >::setup;
      using AmIteratorTmpl<Iterator,DArray<double> >::solve;
      using AmIteratorTmpl<Iterator,DArray<double> >::setClassName;
      
   protected:
  
      // Inherited protected members 
      using Iterator::system;

   private:
      
      /**
      * Find L2 norm of a residual vector.
      */
      double findNorm(DArray<double> const & hist);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double findMaxAbs(DArray<double> const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DArray<double> > & basis, 
                       RingBuffer<DArray<double> > const & hists);

      /**
      * Compute the dot product for one element of the U matrix.
      * 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<DArray<double> > const & resBasis, 
                             int m, int n);

      /**
      * Compute the dot product for one element of the v vector.
      * 
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors
      * \param m row index of the v vector
      */
      double computeVDotProd(DArray<double> const & resCurrent, 
                             RingBuffer<DArray<double> > const & resBasis, 
                             int m);

      /**
      * Update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of past states
      */
      void updateU(DMatrix<double> & U, 
                   RingBuffer<DArray<double> > const & resBasis, 
                   int nHist);

      /**
      * Update the v vector.
      * 
      * \param v v vector
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of past states 
      */
      void updateV(DArray<double> & v, 
                   DArray<double> const & resCurrent, 
                   RingBuffer<DArray<double> > const & resBasis, 
                   int nHist);

      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Add linear combination of basis vectors to trial field.
      * 
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(DArray<double>& trial, 
                        RingBuffer<DArray<double> > const & basis, 
                        DArray<double> coeffs, 
                        int nHist);

      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(DArray<double>& fieldTrial, 
                             DArray<double> const & resTrial, 
                             double lambda);
                             
      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /** 
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements();

      /**
      * Gets the current field vector from the system.
      * 
      * \param curr current field vector
      */ 
      void getCurrent(DArray<double>& curr);

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate();

      /**
      * Compute the residual vector.
      *
      * \param resid current residual vector value
      */
      void getResidual(DArray<double>& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(DArray<double>& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();
      
   
      /**
      * Return true iff all species are treated in closed ensemble.
      */   
      bool isCanonical();

   };

} // namespace Fd1d
} // namespace Pscf
#endif
