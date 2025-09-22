#ifndef R1D_AM_ITERATOR_H
#define R1D_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorDArray.h>  // base class template
#include "Iterator.h"                        // base class template param
#include <pscf/iterator/AmbdInteraction.h>   // member

namespace Pscf {
namespace R1d {

   // Forward declaration
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
   class AmIterator : public AmIteratorDArray<Iterator>
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

   protected:

      // Inherited protected members
      // using Iterator::system;

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

   };

} // namespace R1d
} // namespace Pscf
#endif
