#ifndef RPC_AM_COMPRESSOR_H
#define RPC_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"
#include <prdc/cpu/RField.h>
#include <pscf/iterator/AmIteratorTmpl.h>

namespace Pscf {
namespace Rpc
{

   // Forward declaration
   template <int D> class System;

   // Namespaces from which names can be used without qualification
   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Anderson Mixing compressor.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor
      : public AmIteratorTmpl<Compressor<D>, DArray<double> >
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

      // Inherited public member functions
      using AmIteratorTmpl<Compressor<D>, DArray<double> >::setClassName;

   protected:

      // Inherited protected members
      using ParamComposite::readOptional;
      using Compressor<D>::mdeCounter_;

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

      /**
      * New Basis variable used in updateBasis function
      */
      DArray<double> newBasis_;

      #if 0
      /**
      * Update the basis for residual or field vectors.
      *
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DArray<double> > & basis,
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

      /**
      * Compute mixing parameter lambda
      */
      double computeLambda(double r) override;

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

   };

   // Explicit instantiation declarations
   extern template class AmCompressor<1>;
   extern template class AmCompressor<2>;
   extern template class AmCompressor<3>;

} // namespace Rpc
} // namespace Pscf
#endif
