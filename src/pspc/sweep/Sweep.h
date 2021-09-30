#ifndef PSPC_SWEEP_H
#define PSPC_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/sweep/SweepTmpl.h>        // base class
#include <util/param/ParamComposite.h>        // base class
#include <util/containers/DArray.h>        // base class
#include <util/global.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Solve a sequence of problems along a line in parameter space.
   */
   class Sweep : public SweepTmpl<DArray<double>>
   {

   public:

      /**
      * Default Constructor.
      */
      Sweep();

      /**
      * Destructor.
      */
      ~Sweep();

   protected:

      /**
      * Setup operation at beginning sweep.
      *
      * Must call initializeHistory.
      */
      virtual void setup();

      /**
      * Set non-adjustable system parameters to new values.
      *
      * \param s path length coordinate, in range [0,1]
      */
      virtual void setParameters(double s);

      /**
      * Create guess for adjustable variables by continuation.
      */
      virtual void setGuess(double s);

      /**
      * Call current iterator to solve SCFT problem.
      *
      * Return 0 for sucessful solution, 1 on failure to converge.
      */
      virtual int solve(bool isContinuation);

      /**
      * Reset system to previous solution after iterature failure.
      *
      * The implementation of this function should reset the system state
      * to correspond to that stored in state(0).
      */
      virtual void reset();

      /**
      * Update state(0) and output data after successful convergence
      *
      * The implementation of this function should copy the current 
      * system state into state(0) and output any desired information
      * about the current converged solution.
      */
      virtual void getSolution();

   private:

      DArray<double> stateA_;
      DArray<double> stateB_;
      // DArray<double> stateC_;

   };

} // namespace Pspc
} // namespace Pscf
#endif
