#ifndef RP_AM_COMPRESSOR_H
#define RP_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>   // base class template
#include <util/containers/DArray.h>         // member

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   // Namespaces that can be used implicitly
   using namespace Util;

   /**
   * Anderson mixing compressor.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, both also named AmCompressor,
   * that are defined in Rpc and Rpg namespaces and used in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension of space (D=1, 2, or 3)
   *    - T : Types class (Rpc::Types<D> or Rpg::Types<D>)
   *    - V : 1D array type used for state and residual vectors
   *
   * \see \ref rp_AmCompressor_page "Manual Page"
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T, class V>
   class AmCompressor 
    : public AmIteratorTmpl< typename T::Compressor, V >
   {

   public:

      /// Type for state and residual vectors.
      using VectorT = V;

      /**
      * Read all parameters and initialize.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at
      * the beginning of iteration
      *
      * \param isContinuation true iff continuation within a sweep
      */
      void setup(bool isContinuation) override;

      /**
      * Compress to obtain partial saddle point w+
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress() override;

      /**
      * Return compressor times contributions.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const override;

      /**
      * Clear all timers and mde counter.
      */
      void clearTimers() override;

   protected:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      AmCompressor(System<D,T>& system);

      /**
      * Destructor.
      */
      ~AmCompressor() = default;

      /// Compressor type.
      using CompressorT = typename T::Compressor;

      // Inherited protected member function
      using CompressorT::system;

   private:

      /**
      * Initial values of all fields.
      */
      DArray< typename T::RField > w0_;

      /**
      * Temporary array of w fields used in update function.
      */
      DArray< typename T::RField > wFieldTmp_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

      // Private virtual AM algorithm functions.

      /**
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess() override;

      /**
      * Gets the current field vector from the system.
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
      * Update the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Output relevant system details to the iteration log.
      */
      void outputToLog() override;

      /// Typename alias for base class.
      using AmTmpl = AmIteratorTmpl< CompressorT, VectorT >;

   };

} // namespace Rp
} // namespace Pscf
#endif
