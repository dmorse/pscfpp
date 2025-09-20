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
   * Anderson-Mixing iterator for 1D SCFT.
   *
   * \see \ref r1d_AmIterator_page      "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page "AM Iterator Algorithm"
   *
   * \ingroup R1d_Iterator_Module
   */
   class AmIterator : public AmIteratorTmpl<Iterator, DArray<double> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System 
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
      void readParameters(std::istream& in) override;

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
      void setup(bool isContinuation) override;

   private:

      /// Local copy of interaction, for use with AMBD residual definition
      AmbdInteraction interaction_;

      // -- Virtual functions used to implement AM algorithm -- //

      #if 0
      /**
      * Update the basis for residual or field vectors.
      *
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void 
      updateBasis(RingBuffer<DArray<double> > & basis,
                  RingBuffer<DArray<double> > const & hists) override;

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
                        int nHist) override;

      /**
      * Add predicted error to field trial.
      *
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing
      */
      void addPredictedError(DArray<double>& fieldTrial,
                             DArray<double> const & resTrial,
                             double lambda) override;
      #endif

      // -- Pure virtual functions to exchange data with parent System -- //

      /**
      * Checks if the system has an initial guess.
      */
      bool hasInitialGuess() override;

      /**
      * Compute and returns the number residuals and unknowns.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Gets the current field vector from the system.
      *
      * \param curr current field vector
      */
      void getCurrent(DArray<double>& curr) override;

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate() override;

      /**
      * Compute the residual vector.
      *
      * \param resid current residual vector value
      */
      void getResidual(DArray<double>& resid) override;

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(DArray<double>& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;

      // Non-virtual private function

      /**
      * Return true iff all species are treated in closed ensemble.
      */
      bool isCanonical();

   };

} // namespace R1d
} // namespace Pscf
#endif
