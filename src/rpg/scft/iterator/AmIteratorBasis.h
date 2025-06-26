#ifndef RPG_AM_ITERATOR_BASIS_H
#define RPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"                        // base class argument
#include <pscf/iterator/AmIteratorTmpl.h>    // base class template
#include <pscf/iterator/AmbdInteraction.h>   // member variable
#include <util/containers/DArray.h>          // base class argument
#include <util/containers/RingBuffer.h>      // method input variable

namespace Pscf {
namespace Rpg
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Rpg implementation of the Anderson Mixing iterator.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis
      : public AmIteratorTmpl< Iterator<D>, DArray<double> >
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
      void outputTimers(std::ostream& out);

      // Inherited public member functions
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::solve;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::clearTimers;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;
      using Iterator<D>::nFlexibleParams;

   protected:

      // Inherited protected members
      using ParamComposite::readOptional;
      using ParamComposite::readParamCompositeOptional;
      using ParamComposite::readOptionalFSArray;
      using ParamComposite::setClassName;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::verbose;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::residual;
      using Iterator<D>::system;
      using Iterator<D>::isSymmetric_;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::flexibleParams_;

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation);

   private:

      /// Local copy of interaction, adapted for use AMBD residual definition
      AmbdInteraction interaction_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      /**
      * Set a vector equal to another (assign a = b)
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Compute the inner product of two real vectors.
      */
      double dotProduct(DArray<double> const & a, DArray<double> const & b);

      /**
      * Find the maximum magnitude element of a vector.
      */
      double maxAbs(DArray<double> const & hist);

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of pase residual or field vectors
      */
      void updateBasis(RingBuffer< DArray<double> > & basis, 
                       RingBuffer< DArray<double> > const & hists);

      /**
      * Compute trial field so as to minimize L2 norm of residual.
      * 
      * \param trial resulting trial field (output)
      * \param basis RingBuffer of residual basis vectors.
      * \param coeffs coefficients of basis vectors
      * \param nHist number of prior states stored
      */
      void addHistories(DArray<double>& trial, 
                        RingBuffer<DArray<double> > const & basis, 
                        DArray<double> coeffs, int nHist);

      /**
      * Add predicted error to the trial field.
      * 
      * \param fieldTrial trial field (input/output)
      * \param resTrial predicted error for current trial field
      * \param lambda Anderson-Mixing mixing parameter 
      */
      void addPredictedError(DArray<double>& fieldTrial, 
                             DArray<double> const & resTrial, 
                             double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
     
      /** 
      * Compute the number of elements in the residual vector.
      */
      int nElements();

      /**
      * Get the current w fields and lattice parameters.
      *
      * \param curr current field vector (output)
      */
      void getCurrent(DArray<double>& curr);

      /**
      * Solve MDE for current state of system.
      */
      void evaluate();

      /**
      * Gets the residual vector from system.
      *  
      * \param resid current residual vector (output)
      */
      void getResidual(DArray<double>& resid);

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess trial field configuration
      */
      void update(DArray<double>& newGuess);

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog();

   };

   #ifndef RPG_AM_ITERATOR_BASIS_TPP
   // Suppress implicit instantiation
   extern template class AmIteratorBasis<1>;
   extern template class AmIteratorBasis<2>;
   extern template class AmIteratorBasis<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
#endif
