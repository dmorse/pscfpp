#ifndef PSPC_AM_ITERATOR_H
#define PSPC_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"                     // base class
#include <pscf/iterator/AmIteratorTmpl.h> // base class
#include <util/global.h>                  

namespace Pscf {
namespace Pspc {

   template <int D>
   class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspc_Iterator_Module
   */

   template <int D>
   class AmIterator : public Iterator<D>, public AmIteratorTmpl<DArray<double>>
   {

   public:

      /**
      * Constructor.
      * 
      * \param system system object by reference
      */
      AmIterator(System<D>& system);

      /**
      * Destructor.
      */
      ~AmIterator();

      using AmIteratorTmpl<DArray<double>>::setup;
      using AmIteratorTmpl<DArray<double>>::solve;
      using AmIteratorTmpl<DArray<double>>::readParameters;

   private:

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated
      int nElements();

      /// Gets a reference to the current state of the system
      void getCurrent(DArray<double>& curr);

      /// Runs calculation to evaluate function for fixed point.
      void evaluate();

      /// Gets residual values from system
      void getResidual(DArray<double>& resid);

      /// Updates the system with a passed in state of the iterator.
      void update(DArray<double>& newGuess);

      double findNorm(DArray<double> const & hist);

      double findMaxAbs(DArray<double> const & hist);

      void updateBasis(RingBuffer<DArray<double>> & basis, RingBuffer<DArray<double>> const & hists);

      double computeUDotProd(RingBuffer<DArray<double>> const & resBasis, int m, int n);

      double computeVDotProd(DArray<double> const & resCurrent, RingBuffer<DArray<double>> const & resBasis, int m);

      void updateU(DMatrix<double> & U, RingBuffer<DArray<double>> const & resBasis, int nHist);

      void updateV(DArray<double> & v, DArray<double> const & resCurrent, RingBuffer<DArray<double>> const & resBasis, int nHist);

      void setEqual(DArray<double>& a, DArray<double> const & b);

      void addHistories(DArray<double>& trial, RingBuffer<DArray<double>> const & basis, DArray<double> coeffs, int nHist_);

      void addPredictedError(DArray<double>& fieldTrial, DArray<double> const & resTrial, double lambda);

   };

} // namespace Pspc
} // namespace Pscf
#endif
