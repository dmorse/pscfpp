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

      // Non-virtual private function

      /**
      * Return true iff all species are treated in closed ensemble.
      */
      bool isCanonical();

      // -- Private virtual functions that interact with parent System -- //

      /**
      * Compute and returns the number residuals and unknowns.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Checks if the system has an initial guess.
      */
      bool hasInitialGuess() override;

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

      // --- Private virtual functions for vector math --- //

      /**
      * Assignment for vectors of type T.
      *
      * This function must perform an assignment a = b.
      *
      * \param a  vector to be set (lhs of assignment)
      * \param b  vector value to assign (rhs of assignment)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b) 
      override;

      /**
      * Compute the inner product of two vectors.
      *
      * \param a first vector
      * \param b second vector
      */
      double dotProduct(DArray<double> const & a, DArray<double> const & b) 
      override;

      /**
      * Return the maximum magnitude element of a vector.
      *
      * \param hist  input vector
      */
      double maxAbs(DArray<double> const & hist) 
      override;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a result vector (LHS)
      * \param b first vector (RHS)
      * \param c second vector (RHS)
      */
      void subVV(DArray<double>& a, 
                 DArray<double> const & b, 
		 DArray<double> const & c) 
      override;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a result vector (LHS)
      * \param b input vector (RHS)
      * \param c scalar coefficient (RHS)
      */
      void addEqVc(DArray<double>& a, DArray<double> const & b, double c) 
      override;

   };

} // namespace R1d
} // namespace Pscf
#endif
