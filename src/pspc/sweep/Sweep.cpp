/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include <pscf/sweep/SweepTmpl.tpp>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Setup operation at beginning sweep.
   *
   * Must call initializeHistory.
   */
   void Sweep::setup() 
   {};

   /**
   * Set non-adjustable system parameters to new values.
   *
   * \param s path length coordinate, in range [0,1]
   */
   void Sweep::setParameters(double s) 
   {};

   /**
   * Create guess for adjustable variables by continuation.
   */
   void Sweep::setGuess(double s) 
   {};

   /**
   * Call current iterator to solve SCFT problem.
   *
   * Return 0 for sucessful solution, 1 on failure to converge.
   */
   int solve() 
   {  
      return 0; 
   };

   /**
   * Reset system to previous solution after iterature failure.
   *
   * The implementation of this function should reset the system state
   * to correspond to that stored in state(0).
   */
   void Sweep::reset() 
   {};

   /**
   * Update state(0) and output data after successful convergence
   *
   * The implementation of this function should copy the current 
   * system state into state(0) and output any desired information
   * about the current converged solution.
   */
   void Sweep::getSolution() 
   {};

} // namespace Pspc
} // namespace Pscf
