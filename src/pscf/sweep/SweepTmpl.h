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
   * Solve a sequence of problems along a path through parameter space.
   *
   * \ingroup Pscf_Sweep_Module
   */
   template <typename State>
   class SweepTmpl : public ParamComposite
   {

   public:

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
      * Constructor (protected).
      *
      * The value of historyCapacity depends on the order of continuation,
      * e.g., 2 for 1st order or linear continuation or 3 for 2nd order
      * or quadratic contination. It is passed to this template by a
      * subclass via a protected constructor in order to allow different 
      * derived classes to use different orders of continuation.
      */
      SweepTmpl(int historyCapacity);

      /**
      * Get reference to a stored state, with i=0 being most recent.
      *
      * Call state(i) to return the ith from most recent converged solution.
      */
      State& state(int i)
      {
         UTIL_CHECK(i < historySize_); 
         return *stateHistory_[i]; 
      }

      /**
      * Get the value of s for a stored solution, with i=0 most recent.
      *
      * Call s(i) to return the value of the contour variable s for the 
      * ith from most recent solution.
      */
      double s(int i)
      {
         UTIL_CHECK(i < historySize_); 
         return sHistory_[i]; 
      }

      /**
      * Get the current number of stored previous states.
      */ 
      int historySize()
      {  return historySize_; }

      /**
      * Get the maximum number of stored previous states.
      *
      * The value of historyCapacity is a constant that depends on
      * the order of continuation (i.,e 2 for 1st order continuation,
      * or 3 for 2nd order), and passed as a parameter to the SweepTmpl 
      * constructor.
      */ 
      int historyCapacity()
      {  return historyCapacity_; }

      /**
      * Get the number of converged solutions accepted thus far.
      */ 
      int nAccept()
      {  return nAccept_; }

      /**
      * Initialize variables that track history of solutions.
      *
      * This must be called within the virtual setup() function.
      */
      void initialize();

      /**
      * Check allocation of one state, allocate if necessary.
      *
      * This virtual function is called by SweepTmpl::initialize() 
      * during setup before a sweep to check if all State objects have
      * memory allocated for fields, and allocate if necessary.
      *
      * \param state one stored state of the system.
      */
      virtual void checkAllocation(State & state) = 0;

      /**
      * Setup operation at the beginning of a sweep.
      *
      * Implementations of this function must call initialize().
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
      * The implementation of this function should reset the system 
      * state to correspond to that stored in state(0).
      */
      virtual void reset() = 0;

      /**
      * Update state(0) and output data after successful solution.
      *
      * This function is called by accept(). The implementation of this 
      * function should copy the current system state into state(0) and 
      * output any desired information about the current solution. It
      * may also take any other actions that would normally be taken
      * after acceptance of a converted solution. For example, it may
      * examine the rate of convergence of previous solutions in order 
      * to decide whether to update any parameters used by the SCFT
      * iterator.
      */
      virtual void getSolution() = 0;

      /**
      * Clean up operation at the end of a sweep
      *
      * Default implementation is empty.
      */
      virtual void cleanup();

   private:

      /// Array of State objects, not sequential (work space)
      DArray<State> states_;

      /// Values of s associated with previous solutions
      DArray<double> sHistory_;

      /// Pointers to State objects containing old solutions.
      DArray<State*> stateHistory_;

      /// Maximum number of stored previous states.
      int historyCapacity_;

      /// Current number of stored previous states.
      int historySize_;

      /// Number of converged solutions accepted thus far.
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

      /**
      * Default constructor (private, not implemented to prevent use).
      */
      SweepTmpl();

   };

} // namespace Pscf
#include "SweepTmpl.tpp"
#endif
