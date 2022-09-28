#ifndef PSCF_SWEEP_TMPL_H
#define PSCF_SWEEP_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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

      // Constructor is protected (see below).

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

      /// Whether to write real space concentration field files. 
      bool writeRhoRGrid_;

      /**
      * Constructor (protected).
      *
      * The value of historyCapacity depends on the order of continuation,
      * e.g., 2 for 1st order or linear continuation or 3 for 2nd order
      * or quadratic contination. The value passed to this constructor is
      * a default value that may overridden by a optional parameter in
      * the parameter file format implemented in readParam.
      *
      * \param historyCapacity default maximum number of stored states
      */
      SweepTmpl(int historyCapacity);

      /**
      * Get reference to a stored state, with i=0 being most recent.
      *
      * Call state(i) to return the ith from most recent previous state.
      *
      * \param i history index (i=0 is most recent)
      */
      State& state(int i)
      {
         UTIL_CHECK(i < historySize_); 
         return *stateHistory_[i]; 
      }

      /**
      * Get the value of s for a stored solution, with i = 0 most recent.
      *
      * This function returns the value of the contour variable s for a
      * stored state. Call s(i) to get the value of s for the ith most
      * recent state.
      *
      * \param i history index (i = 0 is most the recent converged state)
      */
      double s(int i) const
      {
         UTIL_CHECK(i < historySize_); 
         return sHistory_[i]; 
      }

      /**
      * Get a coefficient of a previous state in a continuation.
      *
      * An extrapolated trial value for each field or other variables 
      * that describes a state is constructed as a linear superposition 
      * of corresponding values in previous states. Coefficient c(i) is
      * the coefficient of state state(i) in this linear superposition, 
      * where i = 0 denotes the most recent accepted solution and 
      * increasing index i corresponds to increasingly far in the past.
      * Valid values of i are in the range 0 <= i < historySize().
      *
      * The function setCoefficients(double sNew) method computes and 
      * stores values coefficients c(0), ..., c(historySize-1) from 
      * values of sNew (the contour variable of the new state) and 
      * previous values of s. These coefficient values can then be 
      * retrieved by this function. 
      *
      * \param i history index (i=0 is most recent)
      */
      double c(int i) const
      {
         UTIL_CHECK(i >= 0); 
         UTIL_CHECK(i < historySize_); 
         return c_[i]; 
      }

      /**
      * Get the current number of stored previous states.
      */ 
      int historySize() const
      {  return historySize_; }

      /**
      * Get the maximum number of stored previous states.
      *
      * The value of historyCapacity is a constant that is one greater
      * than the maximum order of continuation (e.g., 3 for 2nd order 
      * continuation).  The value is set by passing it as an argument 
      * to the constructor, and is constant after construction.
      */ 
      int historyCapacity() const
      {  return historyCapacity_; }

      /**
      * Get the number of converged solutions accepted thus far.
      *
      * This value is reset to zero by the initialize function, which 
      * must be called by the setup function, and is incremented by one
      * by the accept function.
      */ 
      int nAccept() const
      {  return nAccept_; }

      /**
      * Initialize variables that track history of solutions.
      *
      * This *must* be called within implementation of the setup function.
      */
      void initialize();

      /**
      * Check allocation of one state, allocate if necessary.
      *
      * This virtual function is called by SweepTmpl::initialize() during
      * setup before a sweep to check allocation state and/or allocate 
      * memory for fields in all stored State objects.
      *
      * \param state  an object that represents a state of the system
      */
      virtual void checkAllocation(State & state) = 0;

      /**
      * Setup operation at the beginning of a sweep.
      *
      * The implementations of this function must call initialize().
      */
      virtual void setup() = 0;

      /**
      * Set non-adjustable system parameters to new values.
      *
      * This function should set set values for variables that are treated
      * as input parameters by the SCFT solver, such as block polymer
      * block lengths, chi parameters, species volume fractions or 
      * chemical potentials, etc. The function must modify the values 
      * stored in the parent system to values appropriate to a new value
      * of a contour variable value sNew that is passed as a parameter.
      * The functional dependence of parameters on the contour variable
      * over a range [0,1] is defined by the subclass implementation.
      *
      * \param sNew  new value of path length coordinate, in range [0,1]
      */
      virtual void setParameters(double sNew) = 0;

      /**
      * Create initial guess for the next state by extrapolation.
      *
      * This function should set extrapolated values of the variables that 
      * are modified by the iterative SCFT solver, i.e., values of fields
      * (coefficients of basis functions or values grid points) and unit 
      * cell parameters or domain dimensions for problems involving an 
      * adjustable domain. Values should be extrapolated to a new 
      * contour variable sNew by constructing a linear combination of 
      * corresponding values obtained in previous converged states. 
      * After computing the desired extrapolated values, this function
      * must set these values in the parent system.
      *
      * Extrapolated values of adjustable variables at the new contour
      * variable sNew that is passed as a parameter should be computed
      * for each adjustable variable by constructing a Lagrange polynomial 
      * in s that passes through all stored values, and evaluating the 
      * resulting polynomial at sNew.  This yields an expression for
      * the extrapolated value as a linear combination of stored values
      * with coefficients that depend only the values of sNew and the
      * values of s at previous states. This function should call the 
      * setCoefficients function to compute these coefficients.
      *
      * \param sNew  new value of path length coordinate.
      */
      virtual void extrapolate(double sNew) = 0;

      /**
      * Compute coefficients of previous states for continuation.
      *
      * This function must be called by the implementation of extrapolate.
      *
      * \param sNew  new value of path length coordinate.
      */
      void setCoefficients(double sNew);

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
      * function should copy the current system state into state(0),
      * output any desired information about the current solution,
      * and perform any other operations that should be performed 
      * immediately after acceptance of a converged solution. 
      */
      virtual void getSolution() = 0;

      /**
      * Clean up operation at the end of a sweep
      *
      * Empty default implementation.
      */
      virtual void cleanup();

   private:

      /// Array of State objects, not sequential (work space)
      DArray<State> states_;

      /// Values of s associated with previous solutions
      DArray<double> sHistory_;

      /// Pointers to State objects containing old solutions.
      DArray<State*> stateHistory_;

      /// Coefficients for use during continuation
      DArray<double> c_;

      /// Maximum number of stored previous states.
      int historyCapacity_;

      /// Current number of stored previous states.
      int historySize_;

      /// Number of converged solutions accepted thus far.
      int nAccept_;

      /**
      * Accept a new solution, and update history.
      *
      * This function is called by sweep after a converged solution is 
      * obtained at a new value of the contour variable s. The function
      * calls the pure virtual function getSolution() internally to 
      * copy the system state to the most recent stored solution.
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
