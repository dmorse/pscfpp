#ifndef RPG_AM_ITERATOR_GRID_H
#define RPG_AM_ITERATOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"                        // base class argument
#include <prdc/cuda/types.h>                 // base class argument
#include <pscf/cuda/DeviceArray.h>           // base class argument
#include <pscf/iterator/AmIteratorTmpl.h>    // base class template

#include <pscf/iterator/AmbdInteraction.h>   // member variable
#include <util/containers/DArray.h>          // base class argument
#include <util/containers/RingBuffer.h>      // method input variable

namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

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

      using VectorT = DeviceArray<cudaReal>;
      using Base = AmIteratorTmpl<Iterator<D>, VectorT >;

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
      using Base::solve;
      using Base::clearTimers;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;
      using Iterator<D>::nFlexibleParams;
      using Iterator<D>::stress;

   protected:

      // Inherited protected members
      using Base::verbose;
      using Base::residual;
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

      /// Local copy of interaction, adapted for use AMBD residual definition
      AmbdInteraction interaction_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      // Private virtual functions that interact with parent system

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

      // Private virtual functions for vector math
      
      /**
      * Set vector a equal to vector b (a = b).
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(VectorT& a, VectorT const & b) override;

      /**
      * Compute and return inner product of two real fields.
      */
      double dotProduct(VectorT const & a, 
                        VectorT const & b) override;

      /**
      * Find the maximum magnitude element of a residual vector.
      *  
      * \param a input vector
      */
      double maxAbs(VectorT const & a) override;

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
		   VectorT const & b, double c) override;

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
