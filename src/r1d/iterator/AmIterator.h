#ifndef R1D_AM_ITERATOR_H
#define R1D_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <r1d/solvers/Mixture.h>
#include <pscf/iterator/AmbdInteraction.h>
#include <pscf/iterator/AmIteratorTmpl.h>

namespace Pscf {
namespace R1d
{

   class System;

   using namespace Util;

   /**
   * Anderson-Mixing iterator.
   *
   * \ingroup R1d_Iterator_Module
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
      using AmIteratorTmpl<Iterator,DArray<double> >::solve;
      using AmIteratorTmpl<Iterator,DArray<double> >::setClassName;

   protected:

      // Inherited protected members
      using Iterator::system;

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation);

   private:

      // Local copy of interaction, adapted for use AMBD residual definition
      AmbdInteraction interaction_;

      // -- Virtual functions used to implement AM algorithm -- //

      /**
      * Assign one field to another.
      *
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Find L2 norm of a residual vector.
      */
      double
      dotProduct(DArray<double> const & a, DArray<double> const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(DArray<double> const & hist);

      /**
      * Update the basis for residual or field vectors.
      *
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DArray<double> > & basis,
                       RingBuffer<DArray<double> > const & hists);

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

      // -- Virtual funtions to exchange data with parent System -- //

      /**
      * Checks if the system has an initial guess.
      */
      bool hasInitialGuess();

      /**
      * Compute and returns the number residuals and unknowns.
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

} // namespace R1d
} // namespace Pscf
#endif
