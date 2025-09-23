#ifndef RPG_AM_ITERATOR_BASIS_H
#define RPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"                        // base class argument
#include <pscf/iterator/AmIteratorDArray.h>  // base class template
#include <pscf/iterator/AmbdInteraction.h>   // member variable
#include <util/containers/DArray.h>          // function argument
#include <util/containers/RingBuffer.h>      // function argument

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   //using namespace Prdc;

   /**
   * Anderson Mixing iterator with imposed space-group symmetry.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis : public AmIteratorDArray< Iterator<D> >
   {

   public:

      /**
      * Constructor.
      *   
      * \param system parent system object
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
      void readParameters(std::istream& in);

      /**
      * Output timing results to log file.
      *
      * \param out  output stream for timer report
      */
      void outputTimers(std::ostream& out) const;

      // Inherited public member functions
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::solve;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::clearTimers;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;
      using Iterator<D>::nFlexibleParams;
      using Iterator<D>::stress;

   protected:

      // Inherited protected members
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::verbose;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::residual;
      using Iterator<D>::system;
      using Iterator<D>::isSymmetric_;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::flexibleParams_;
      using ParamComposite::readOptional;
      using ParamComposite::readParamCompositeOptional;
      using ParamComposite::readOptionalFSArray;
      using ParamComposite::setClassName;

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation);

   private:

      /// Local copy of interaction, for use with AMBD residual definition
      AmbdInteraction interaction_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      // Private virtual functions that interact with parent System.
      
      /** 
      * Compute the number of elements in the residual vector.
      */
      int nElements() override;

      /**
      * Check if the system has an initial guess.
      */
      bool hasInitialGuess() override;
     
      /**
      * Get the current w fields and lattice parameters.
      *
      * \param curr current field vector (output)
      */
      void getCurrent(DArray<double>& curr) override;

      /**
      * Solve MDE for current state of system.
      */
      void evaluate() override;

      /**
      * Gets the residual vector from system.
      *  
      * \param resid current residual vector (output)
      */
      void getResidual(DArray<double>& resid) override;

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess trial field configuration
      */
      void update(DArray<double>& newGuess) override;

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog() override;

      #if 0
      // Private virtual functions for vector math
      
      /**
      * Set a vector equal to another (assign a = b)
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(DArray<double>& a, 
                    DArray<double> const & b) override;

      /**
      * Compute the inner product of two real vectors.
      */
      double dotProduct(DArray<double> const & a, 
                        DArray<double> const & b) override;

      /**
      * Find the maximum magnitude element of a vector.
      */
      double maxAbs(DArray<double> const & hist) override;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a result vector (LHS)
      * \param b first vector (RHS)
      * \param c second vector (RHS)
      */
      void subVV(DArray<double>& a, 
                 DArray<double> const & b, 
		 DArray<double> const & c) override;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a result vector (LHS)
      * \param b input vector (RHS)
      * \param c scalar coefficient (RHS)
      */
      void addEqVc(DArray<double>& a, 
		   DArray<double> const & b, double c) override;
      #endif


   };

   // Explicit instantiation declarations
   extern template class AmIteratorBasis<1>;
   extern template class AmIteratorBasis<2>;
   extern template class AmIteratorBasis<3>;

} // namespace Rpg
} // namespace Pscf
#endif
