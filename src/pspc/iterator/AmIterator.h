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

      // Inherited public member functions
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::setup;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::solve;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::setClassName;
      using Iterator<D>::maskField;
      using Iterator<D>::externalField;
      using Iterator<D>::hasMask;
      using Iterator<D>::hasExternalField;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;

   protected:
  
      // Inherited protected members 
      using ParamComposite::readOptional;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::flexibleParams_;

   private:

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;
      
      /**
      * Find L2 norm of a residual vector.
      */
      double findNorm(FieldCPU const & hist);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double findMaxAbs(FieldCPU const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<FieldCPU> & basis, 
                       RingBuffer<FieldCPU> const & hists);

      /**
      * Compute the dot product for one element of the U matrix.
      * 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<FieldCPU> const & resBasis, 
                             int m, int n);

      /**
      * Compute the dot product for one element of the v vector.
      * 
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors
      * \param m row index of the v vector
      */
      double computeVDotProd(FieldCPU const & resCurrent, 
                             RingBuffer<FieldCPU> const & resBasis, 
                             int m);

      /**
      * Update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of past states
      */
      void updateU(DMatrix<double> & U, 
                   RingBuffer<FieldCPU> const & resBasis, 
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
                   FieldCPU const & resCurrent, 
                   RingBuffer<FieldCPU> const & resBasis, 
                   int nHist);

      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(FieldCPU& a, FieldCPU const & b);

      /**
      * Add linear combination of basis vectors to trial field.
      * 
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(FieldCPU& trial, 
                        RingBuffer<FieldCPU> const & basis, 
                        DArray<double> coeffs, 
                        int nHist);

      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(FieldCPU& fieldTrial, 
                             FieldCPU const & resTrial, 
                             double lambda);

      /**
      * Does the system has an initial guess for the field?
      */
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
      void getCurrent(FieldCPU& curr);

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
      void getResidual(FieldCPU& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(FieldCPU& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();

   };

} // namespace Pspc
} // namespace Pscf
#endif
