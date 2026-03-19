#ifndef RPC_AM_ITERATOR_GRID_H
#define RPC_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBase.h"                  // base class
#include <pscf/iterator/AmbdInteraction.h>   // member

namespace Pscf {
namespace Rpc {

   template <int D> class System;

   using namespace Util;

   /**
   * Anderson Mixing iterator on grid (no space-group symmetry).
   *
   * \see \ref rp_AmIteratorGrid_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   * \ingroup Rpc_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorGrid
    : public AmIteratorTmpl< Iterator<D>, DRArray<double> >
   {

   public:

      /// Alias for type of state and residual vectors.
      using VectorT = DRArray<double>;

      /// Alias for base class.
      using AmIterTmplT = AmIteratorTmpl< Iterator<D>, VectorT >;

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      AmIteratorGrid(System<D>& system);

      /**
      * Destructor.
      */
      ~AmIteratorGrid();

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

      // Inherited protected members
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::flexibleParams_;

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation) override;

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
      template <typename T> using HostArrayT = DArray<T>;

   };

   // Explicit instantiation declarations
   extern template class AmIteratorGrid<1>;
   extern template class AmIteratorGrid<2>;
   extern template class AmIteratorGrid<3>;

} // namespace Rpc
} // namespace Pscf
#endif
