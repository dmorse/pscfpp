#ifndef PSPC_AM_ITERATOR_H
#define PSPC_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspc implementation of the Anderson Mixing iterator.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class AmIterator : public AmIteratorTmpl<Iterator<D>, FieldCPU>
   {

   public:
      /**
      * Constructor.
      */
      AmIterator(System<D>& system);

      ~AmIterator();

      using AmIteratorTmpl<Iterator<D>,FieldCPU>::setup;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::solve;
      using AmIteratorTmpl<Iterator<D>,FieldCPU>::readParameters;
      using Iterator<D>::sys_;

   private:

      double findNorm(FieldCPU const & hist);

      double findRelNorm(FieldCPU const & hist);

      double findMaxAbs(FieldCPU const & hist);

      void updateBasis(RingBuffer<FieldCPU> & basis, RingBuffer<FieldCPU> const & hists);

      double computeUDotProd(RingBuffer<FieldCPU> const & resBasis, int m, int n);

      double computeVDotProd(FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int m);

      void updateU(DMatrix<double> & U, RingBuffer<FieldCPU> const & resBasis, int nHist);

      void updateV(DArray<double> & v, FieldCPU const & resCurrent, RingBuffer<FieldCPU> const & resBasis, int nHist);

      void setEqual(FieldCPU& a, FieldCPU const & b);

      void addHistories(FieldCPU& trial, RingBuffer<FieldCPU> const & basis, DArray<double> coeffs, int nHist);

      void addPredictedError(FieldCPU& fieldTrial, FieldCPU const & resTrial, double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated
      int nElements();

      /// Gets a reference to the current state of the system
      void getCurrent(FieldCPU& curr);

      /// Runs calculation to evaluate function for fixed point.
      void evaluate();

      /// Gets residual values from system
      void getResidual(FieldCPU& resid);

      /// Updates the system with a passed in state of the iterator.
      void update(FieldCPU& newGuess);

      /// Outputs relevant system details to the iteration log
      void outputToLog();

   };

} // namespace Pspc
} // namespace Pscf
#endif
