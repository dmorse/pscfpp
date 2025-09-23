#ifndef RPG_AM_COMPRESSOR_H
#define RPG_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                    // base class argument
#include <prdc/cuda/types.h>               // base class argument
#include <pscf/cuda/DeviceArray.h>         // base class argument
#include <pscf/iterator/AmIteratorTmpl.h>  // base class template

#include <prdc/cuda/RField.h>              // member
#include <util/containers/DArray.h>        // member

#include <iostream>

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Anderson Mixing compressor.
   *
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor 
     : public AmIteratorTmpl< Compressor<D>, DeviceArray<cudaReal> >
   {

   public:

      /// Type for state and residual vectors
      using VectorT = DeviceArray<cudaReal>;

      /// Alias for base class
      using Base =  AmIteratorTmpl< Compressor<D>, VectorT >;

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
      void readParameters(std::istream& in);

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at the
      * beginning of iteration
      *
      * \param isContinuation true iff continuation within a sweep
      */
      void setup(bool isContinuation);

      /**
      * Compress to obtain partial saddle point w+
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress();

      /**
      * Compute mixing parameter lambda
      */
      double computeLambda(double r);

      /**
      * Return how many times MDE has been solved.
      */
      int mdeCounter();

      /**
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out) const;

      /**
      * Clear all timers (reset accumulated time to zero).
      */
      void clearTimers();

   protected:

      // Inherited protected members
      using Compressor<D>::mdeCounter_;
      using Compressor<D>::system;
      using ParamComposite::readOptional;

   private:

      /**
      * Count how many times MDE has been solved.
      */
      int counter_;

      /**
      * Current values of the fields
      */
      DArray< RField<D> > w0_;

      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;

      /**
      * Template w Field used in update function
      */
      DArray< RField<D> > wFieldTmp_;

      // Private virtual functions that interact with parent System

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
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;

      // Private virtual functions for vector math

      /**
      * Assign one field to another.
      *
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(VectorT& a,
                    VectorT const & b) override;

      /**
      * Compute the inner product of two vectors
      */
      double dotProduct(VectorT const & a,
                        VectorT const & b) override;

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(VectorT const & hist) override;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a result vector (LHS)
      * \param b first vector (RHS)
      * \param c second vector (RHS)
      */
      void subVV(VectorT& a, 
                 VectorT const & b, 
		 VectorT const & c) override;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a result vector (LHS)
      * \param b input vector (RHS)
      * \param c scalar coefficient (RHS)
      */
      void addEqVc(VectorT& a, 
		   VectorT const & b, 
		   double c) override;

   };

   // Explicit instantiation declarations
   extern template class AmCompressor<1>;
   extern template class AmCompressor<2>;
   extern template class AmCompressor<3>;

} // namespace Rpg
} // namespace Pscf
#endif
