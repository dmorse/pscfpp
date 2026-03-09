#ifndef RPG_AM_ITERATOR_GRID_H
#define RPG_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorDev.h"                   // base class
#include <pscf/iterator/AmbdInteraction.h>   // member variable
#include <util/containers/RingBuffer.h>      // method input variable

namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;

   /**
   * Rpg implementation of the Anderson Mixing iterator.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class AmIteratorGrid 
     : public AmIteratorTmpl<Iterator<D>, DeviceArray<cudaReal> >
   {

   public:

      // Public typename aliases.

      /// Typename for type of state and residual vectors.
      using VectorT = DeviceArray<cudaReal>;

      /// Alias for base class.
      using AmIterTmplT = AmIteratorTmpl<Iterator<D>, VectorT >;

      // Public member functions

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
      void outputTimers(std::ostream& out) const;

      // Inherited public member functions
      using AmIterTmplT::solve;
      using AmIterTmplT::clearTimers;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;
      using Iterator<D>::nFlexibleParams;
      using Iterator<D>::stress;

   protected:

      // Inherited protected members
      using AmIterTmplT::verbose;
      using AmIterTmplT::residual;
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
      void setup(bool isContinuation) override;

   private:

      /// Local copy of interaction, adapted for AMBD residual definition
      AmbdInteraction interaction_;

      /// Scale factor for stress residual elements
      double scaleStress_;

      // Private virtual functions that interact with parent system

      /** 
      * Compute and return the number of elements in the residual vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system have an initial guess for state vector?
      */
      bool hasInitialGuess() override;
     
      /**
      * Get the current state vector (w fields and lattice parameters).
      *
      * \param curr current state vector (output)
      */
      void getCurrent(VectorT& curr) override;

      /**
      * Solve MDE for current state of system.
      */
      void evaluate() override;

      /**
      * Gets the residual vector from system.
      *  
      * \param resid  current residual vector (output)
      */
      void getResidual(VectorT& resid) override;

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess  trial field configuration (output)
      */
      void update(VectorT& newGuess) override;

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog() override;

      // --- Private member function specific to this implementation --- 
      
      /**
      * Calculate the average value of an array.
      * 
      * \param field  input array
      */
      cudaReal findAverage(VectorT const & field);

   };

   // Explicit instantiation declarations
   extern template class AmIteratorGrid<1>;
   extern template class AmIteratorGrid<2>;
   extern template class AmIteratorGrid<3>;

} // namespace Rpg
} // namespace Pscf
#endif
