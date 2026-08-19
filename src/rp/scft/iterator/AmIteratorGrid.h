#ifndef RP_AM_ITERATOR_GRID_H
#define RP_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>   // base class template
#include <rp/scft/iterator/Iterator.h>      // base class argument
#include <pscf/iterator/AmbdInteraction.h>  // member
#include <pscf/backends/TmplDeclare.h>      // declaration macros

#include <iostream>

// Forward declaration
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class RField;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Anderson mixing iterator on grid (no space-group symmetry).
   *
   * This variant of the Anderson mixing algorithm uses a regular 
   * computational mesh to represent all fields, with no imposed symmetry.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_AmIteratorGrid_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   * \ingroup Rp_Scft_Iterator_Module
   */
   template <int D, class T>
   class AmIteratorGrid
    : public AmIteratorTmpl< Iterator<D,T>, typename T::RDevArray >
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      AmIteratorGrid(System<D,T>& system);

      /**
      * Destructor.
      */
      virtual ~AmIteratorGrid();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in) override;

      /**
      * Output timing results to log file.
      *
      * \param out  output stream for timer report
      */
      void outputTimers(std::ostream& out) const override;

   protected:

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation) override;

      /// Alias for type of state and residual vectors.
      using VectorT = typename T::RDevArray;

      /// Alias for base class.
      using AmIterTmplT = AmIteratorTmpl< Iterator<D,T> , VectorT >;

      // Inherited protected members
      using Iterator<D,T>::system;
      using Iterator<D,T>::flexibleParams_;

   private:

      /// Local copy of interaction, adapted for AMBD residual definition
      AmbdInteraction interaction_;

      /// Scale factor for stress residual elements
      double scaleStress_;

      // Private overridden virtual functions

      /**
      * Compute and return the number of elements in the residual vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system have an initial guess for the state vector?
      */
      bool hasInitialGuess() override;

      /**
      * Get the current state vector (w fields and lattice parameters).
      *
      * \param curr  current state vector (output)
      */
      void getCurrent(VectorT& curr) override;

      /**
      * Have the system perform a computation using new state.
      *
      * Solves the modified diffusion equation, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate() override;

      /**
      * Compute the residual vector.
      *
      * \param resid  current residual vector (output)
      */
      void getResidual(VectorT& resid) override;

      /**
      * Update the system with a new state vector.
      *
      * \param newGuess  new state vector (input)
      */
      void update(VectorT& newGuess) override;

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog() override;

      // Private type aliases
      using RealT = double;
      using RFieldT = RField<D,T>;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(AmIteratorGrid)

} // namespace Rp
} // namespace Pscf
#endif
