/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include <pspc/System.h>
#include <pscf/sweep/SweepTmpl.tpp>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Setup operation at beginning sweep.
   *
   * Must call initializeHistory.
   */
   template <int D>
   void Sweep<D>::setup() 
   {
      for (int i=0; i < 2; ++i) {
         state(i).setSystem(system());
      }
      initializeHistory(state(0), state(1));

      trial_.setSystem(system());
   };

   /**
   * Set non-adjustable system parameters to new values.
   *
   * \param s path length coordinate, in range [0,1]
   */
   template <int D>
   void Sweep<D>::setParameters(double s) 
   {
      // Empty default implementation to allow Sweep<D> to be compiled.
   };

   /**
   * Create guess for adjustable variables by continuation.
   */
   template <int D>
   void Sweep<D>::setGuess(double sNew) 
   {
      DArray<double> coeffs_;
      coeffs_.allocate(nHistory);
      // Set coeffs for appropriate continuation order
      // trial_.linearCombination(history(), coeffs_);
      // Set unit cell and wFields in system
   };

   /**
   * Call current iterator to solve SCFT problem.
   *
   * Return 0 for sucessful solution, 1 on failure to converge.
   */
   template <int D>
   int Sweep<D>::solve(bool isContinuation)
   {  return system().iterate(); };

   /**
   * Reset system to previous solution after iterature failure.
   *
   * The implementation of this function should reset the system state
   * to correspond to that stored in state(0).
   */
   template <int D>
   void Sweep<D>::reset()
   {  state(0).setSystemState(); }

   /**
   * Update state(0) and output data after successful convergence
   *
   * The implementation of this function should copy the current 
   * system state into state(0) and output any desired information
   * about the current converged solution.
   */
   template <int D>
   void Sweep<D>::getSolution() 
   { 
      state(0).getSystemState(); 
      // Output operations (see example in Fd1d)
   };

} // namespace Pspc
} // namespace Pscf
