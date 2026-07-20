#ifndef RP_AM_ITERATOR_BASIS_H
#define RP_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>    // base class template
#include <rp/scft/iterator/Iterator.h>       // base template argument
#include <util/containers/DArray.h>          // base template argument
#include <pscf/iterator/AmbdInteraction.h>   // direct base member 

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Anderson Mixing iterator with imposed space-group symmetry.
   * 
   * This variant of the Anderson mixing algorithm uses an expansion in
   * symmetry-adapted Fourier basis functions to represent all fields,
   * thus automatically imposing a user-specified space group symmetry.
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named AmIteratorBasis,
   * that are defined in Rpc and Rpg namespaces and used in the pscf_rpc
   * and pscf_rpg programs, respectively. 
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Cpp<D> or CudaTp<D>
   *
   * \see \ref rp_AmIteratorBasis_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   *
   * \ingroup Rp_Scft_Iterator_Module
   */
   template <int D, class T>
   class AmIteratorBasis 
    : public AmIteratorTmpl<Iterator<D,T>, DArray<double> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system System object associated with this iterator.
      */
      AmIteratorBasis(System<D,T>& system);

      /**
      * Destructor.
      */
      virtual ~AmIteratorBasis();

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

      /// Alias for type of residual and state vectors
      using VectorT = DArray<double>;

      /// Alias for base class.
      using AmIterTmplT = AmIteratorTmpl< Iterator<D,T>, VectorT>;

      // Inherited protected members
      using Iterator<D,T>::system;
      using Iterator<D,T>::flexibleParams_;

   private:

      /// Local copy of interaction, adapted for AMBD residual definition
      AmbdInteraction interaction_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      // Private virtual functions that interact with parent system

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess() override;

      /**
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Get the current w fields and lattice parameters.
      *
      * \param curr current field vector
      */
      void getCurrent(VectorT& curr) override;

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
      void getResidual(VectorT& resid) override;

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;

   };

} // namespace Rp
} // namespace Pscf
#endif
