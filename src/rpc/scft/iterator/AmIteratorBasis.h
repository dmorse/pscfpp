#ifndef RPC_AM_ITERATOR_BASIS_H
#define RPC_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"                            // base class argument
#include <pscf/iterator/AmIteratorTmpl.h>        // base class template
#include <pscf/iterator/AmbdInteraction.h>       // member variable
#include <util/containers/DArray.h>              // base class argument
#include <util/containers/RingBuffer.h>          // method input variable

namespace Pscf {
namespace Rpc
{

   template <int D> class System;

   using namespace Util;

   /**
   * Rpc implementation of the Anderson Mixing iterator with symmetry.
   * 
   * \see \ref rpc_AmIteratorBasis_page "Manual Page"
   * \see \ref pscf_AmIteratorTmpl_page  "AM Iteration Algorithm"
   *
   * \ingroup Rpc_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis
      : public AmIteratorTmpl< Iterator<D>, DArray<double> >
   {

   public:

      using Base = AmIteratorTmpl< Iterator<D>, DArray<double> >;

      /**
      * Constructor.
      *
      * \param system System object associated with this iterator.
      */
      AmIteratorBasis(System<D>& system);

      /**
      * Destructor.
      */
      ~AmIteratorBasis();

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

      // Inherited public member functions
      using Base::solve;
      using Base::clearTimers;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;
      using Iterator<D>::nFlexibleParams;
      using Iterator<D>::stress;

   protected:

      // Inherited protected members
      using ParamComposite::setClassName;
      using ParamComposite::readOptional;
      using ParamComposite::readParamCompositeOptional;
      using ParamComposite::readOptionalFSArray;
      using Base::verbose;
      using Base::residual;
      using Iterator<D>::system;
      using Iterator<D>::isSymmetric_;
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

      // Private virtual functions for vector math

      /**
      * Vector assignment, a = b.
      *
      * This function must perform an assignment a = b.
      *
      * \param a  vector to be set (LHS)
      * \param b  vector value to assign (RHS)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b) 
      override;

      /**
      * Compute and return the inner product of two vectors.
      *
      * \param a first vector
      * \param b second vector
      */
      double dotProduct(DArray<double> const & a, 
                        DArray<double> const & b) 
      override;

      /**
      * Return the maximum magnitude element of a vector.
      *
      * \param hist  input vector
      */
      virtual double maxAbs(DArray<double> const & hist) override;

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
      void addEqVc(DArray<double>& a, 
		   DArray<double> const & b, double c) 
      override;

   };

   // Explicit instantiation declarations
   extern template class AmIteratorBasis<1>;
   extern template class AmIteratorBasis<2>;
   extern template class AmIteratorBasis<3>;

} // namespace Rpc
} // namespace Pscf
#endif
