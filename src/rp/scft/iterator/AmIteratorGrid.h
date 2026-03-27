#ifndef RP_AM_ITERATOR_GRID_H
#define RP_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>    // base class
#include <pscf/iterator/AmbdInteraction.h>   // member
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Anderson mixing iterator on grid (no space-group symmetry).
   *
   * This variant of the Anderson mixing algorithm uses a regular 
   * computational mesh to represent all fields, with no imposed symmetry.
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named AmIteratorGrid,
   * that are defined in Rpc and Rpg namespaces and used in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * Typename T::Vector must be Util::DArray<double> for use in the
   * Rpc namespace and Pscf::DevArray<cudaReal> for use in the Rpg
   * namespace. Both classes allow the same container to either own
   * data or be associated with a slice of data that is owned by 
   * another instance of the same class.
   *
   * \see \ref rp_AmIteratorGrid_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   * \ingroup Rp_Scft_Iterator_Module
   */
   template <int D, class T>
   class AmIteratorGrid
    : public AmIteratorTmpl< typename T::Iterator, typename T::Vector >
   {

   public:

      /// Alias for Iterator type.
      using IteratorT = typename T::Iterator;

      /// Alias for type of state and residual vectors.
      using VectorT = typename T::Vector;

      /// Alias for base class.
      using AmIterTmplT = AmIteratorTmpl< IteratorT, VectorT >;

      /**
      * Constructor.
      *
      * \param system  parent system object
      */
      AmIteratorGrid(typename T::System& system);

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
      using IteratorT::system;
      using IteratorT::flexibleParams_;

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
      using RFieldT = typename T::RField;

      //template <typename Data> 
      //using HostArrayT = (typename T::HostDArray)<Data>;

   };

} // namespace Rp
} // namespace Pscf
#endif
