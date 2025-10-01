#ifndef RPC_AM_COMPRESSOR_H
#define RPC_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                     // base class template argument
#include <prdc/cpu/RField.h>                // member
#include <pscf/iterator/AmIteratorDArray.h> // base class template
#include <util/containers/DArray.h>         // member

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class System;

   // Namespaces from which names can be used without qualification
   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Anderson mixing compressor.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor
      : public AmIteratorDArray< Compressor<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system System object associated with this compressor.
      */
      AmCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~AmCompressor();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in) override;

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at the
      * beginning of iteration
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
      */
      void outputTimers(std::ostream& out) const override;

      /**
      * Clear all timers (reset accumulated time to zero).
      */
      void clearTimers() override;

   private:

      /**
      * How many times MDE has been solved for this stochastic move?
      */
      int itr_;

      /**
      * Current values of the fields
      */
      DArray< RField<D> > w0_;

      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;

      /**
      * Temporary w field used in update function
      */
      DArray< RField<D> > wFieldTmp_;

      // Private virtual functions that interact with parent system

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

      // Inherited private members
      using Compressor<D>::system;

      // Indirect (grandparent) base class.
      using AmTmpl = AmIteratorTmpl< Compressor<D>, DArray<double> >;

   };

   // Explicit instantiation declarations
   extern template class AmCompressor<1>;
   extern template class AmCompressor<2>;
   extern template class AmCompressor<3>;

} // namespace Rpc
} // namespace Pscf
#endif
