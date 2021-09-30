#ifndef PSCF_SWEEP_TMPL_H
#define PSCF_SWEEP_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>                     // base class

namespace Pscf {

   using namespace Util;

   /**
   * Solve a sequence of problems along a line in parameter space.
   *
   * \ingroup Pscf_Sweep_Module
   */
   template <typename State>
   class SweepTmpl : public ParamComposite
   {

   static const int nHistory = 2;

   public:

      /**
      * Default Constructor.
      */
      SweepTmpl();

      /**
      * Destructor.
      */
      ~SweepTmpl();

      /**
      * Read ns and baseFileName parameters.
      *
      * \param in input stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Iterate to solution.
      */
      virtual void sweep();

   protected:

      /// Number of steps. 
      int ns_;

      /// Base name for output files
      std::string baseFileName_;

      /**
      * Get reference to stored solution, with i=0 being most recent.
      *
      * Call state(i) to return the ith from most recent converged solution.
      */
      State& state(int i)
      { return *stateHistory_[i]; }

      /**
      * Get the value of s for a stored solution, with i=0 most recent.
      *
      * Call s(i) to return the ith from most recent solution.
      */
      double s(int i)
      { return sHistory_[i]; }

      /**
      * Get the number of stored solutions.
      */ 
      int historySize()
      {  return historySize_; }

      /**
      * Get the number of converged solutions accepted thus far in this sweep.
      */ 
      int nAccept()
      {  return nAccept_; }

      /**
      * Initialize variables that track history of solutions.
      *
      * This must be called within the setup() function.
      */
      void initializeHistory(State& state0, State& state1);

      /**
      * Setup operation at the beginning of a sweep.
      *
      * Implementations of this function must call initializeHistory.
      */
      virtual void setup() = 0;

      /**
      * Set non-adjustable system parameters to new values.
      *
      * \param s path length coordinate, in range [0,1]
      */
      virtual void setParameters(double s) = 0;

      /**
      * Create guess for adjustable variables by continuation.
      */
      virtual void setGuess(double sNew) = 0;

      /**
      * Call current iterator to solve SCFT problem.
      *
      * Return 0 for sucessful solution, 1 on failure to converge.
      */
      virtual int solve(bool isContinuation) = 0;

      /**
      * Reset system to previous solution after iterature failure.
      *
      * The implementation of this function should reset the system state
      * to correspond to that stored in state(0).
      */
      virtual void reset() = 0;

      /**
      * Update state(0) and output data after successful solution.
      *
      * This function is called by accept(). The implementation of this 
      * function should copy the current system state into state(0) and 
      * output any desired information about the current solution.
      */
      virtual void getSolution() = 0;

   private:

      // Values of s associated with previous solutions
      FArray<double, nHistory> sHistory_;

      /// Pointers State objects containing old solutions.
      FArray<State*, nHistory> stateHistory_;

      // Number of previous solutions currently stored.
      int historySize_;

      // Number of converged solutions accepted thus far.
      int nAccept_;

      /**
      * Accept a new solution, and update history.
      *
      * This function should be called whenever a converted solution is 
      * obtained. The solution calls getSolution() internally to copy
      * the system state to the most recent stored solution.
      *
      * \param s value of s associated with new solution.
      */
      void accept(double s);

   };

} // namespace Pscf
#include "SweepTmpl.tpp"
#endif
