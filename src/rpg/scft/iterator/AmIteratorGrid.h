#ifndef RPG_AM_ITERATOR_GRID_H
#define RPG_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"                        // base class
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
   class AmIteratorGrid : public AmIteratorTmpl<Iterator<D>, FieldCUDA>
   {

   public:

      /**
      * Constructor.
      *   
      * \param system parent system object
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
      void readParameters(std::istream& in);

      /**
      * Output timing results to log file.
      *
      * \param out  output stream for timer report
      */
      void outputTimers(std::ostream& out);

      // Inherited public member functions
      using AmIteratorTmpl<Iterator<D>, FieldCUDA>::solve;
      using AmIteratorTmpl<Iterator<D>, FieldCUDA>::clearTimers;
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
      using AmIteratorTmpl<Iterator<D>, FieldCUDA>::verbose;
      using AmIteratorTmpl<Iterator<D>, FieldCUDA>::residual;
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

      // Virtual functions used to implement AM algorithm
      
      /**
      * Set vector a equal to vector b (a = b).
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(FieldCUDA& a, FieldCUDA const & b);

      /**
      * Compute and return inner product of two real fields.
      */
      double dotProduct(FieldCUDA const & a, FieldCUDA const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      *  
      * \param a input vector
      */
      double maxAbs(FieldCUDA const & a);

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer of basis vectors.
      * \param hists RingBuffer of previous vectors.
      */
      void updateBasis(RingBuffer<FieldCUDA> & basis, 
                       RingBuffer<FieldCUDA> const & hists);

      /**
      * Compute trial field so as to minimize L2 norm of residual.
      * 
      * \param trial resulting trial field (output)
      * \param basis RingBuffer of residual basis vectors.
      * \param coeffs coefficients of basis vectors
      * \param nHist number of prior states stored
      */
      void addHistories(FieldCUDA& trial, 
                        RingBuffer<FieldCUDA> const & basis, 
                        DArray<double> coeffs, int nHist);

      /**
      * Add predicted error to the trial field.
      * 
      * \param fieldTrial trial field (input/output)
      * \param resTrial predicted error for current trial field
      * \param lambda Anderson-Mixing mixing parameter 
      */
      void addPredictedError(FieldCUDA& fieldTrial, 
                             FieldCUDA const & resTrial, 
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
      void getCurrent(FieldCUDA& curr);

      /**
      * Solve MDE for current state of system.
      */
      void evaluate();

      /**
      * Gets the residual vector from system.
      *  
      * \param resid current residual vector (output)
      */
      void getResidual(FieldCUDA& resid);

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess trial field configuration
      */
      void update(FieldCUDA& newGuess);

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog();

      // --- Private member functions specific to this implementation --- 
      
      /**
      * Calculate the average value of an array.
      * 
      * \param field  input array
      */
      cudaReal findAverage(FieldCUDA const & field);

   };

   #ifndef RPG_AM_ITERATOR_GRID_TPP
   // Suppress implicit instantiation
   extern template class AmIteratorGrid<1>;
   extern template class AmIteratorGrid<2>;
   extern template class AmIteratorGrid<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
#endif
