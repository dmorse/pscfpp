#ifndef PSCF_SWEEP_TMPL_TPP
#define PSCF_SWEEP_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepTmpl.h"

namespace Pscf {

   using namespace Util;

   /*
   * Constructor
   */
   template <class State>
   SweepTmpl<State>::SweepTmpl(int historyCapacity)
    : ns_(0),
      baseFileName_(),
      historyCapacity_(historyCapacity)
   {
      setClassName("SweepTmpl"); 
      states_.allocate(historyCapacity_);
      stateHistory_.allocate(historyCapacity_);
      sHistory_.allocate(historyCapacity_);
      c_.allocate(historyCapacity_);
   }

   /*
   * Destructor.
   */
   template <class State>
   SweepTmpl<State>::~SweepTmpl()
   {}

   /*
   * Read parameters.
   */
   template <class State>
   void SweepTmpl<State>::readParameters(std::istream& in)
   {
      read<int>(in, "ns", ns_);
      read<std::string>(in, "baseFileName", baseFileName_);
   }

   template <class State>
   void SweepTmpl<State>::sweep()
   {

      // Compute and output ds
      double ds = 1.0/double(ns_);
      double ds0 = ds;
      std::cout << std::endl;
      std::cout << "ns = " << ns_ << std::endl;
      std::cout << "ds = " << ds  << std::endl;

      // Initial setup, before sweep
      setup();

      // Solve for initial state of sweep
      int error;
      double sNew = 0.0;
      std::cout << std::endl;
      std::cout << "Attempt s = " << sNew << std::endl;
      bool isContinuation = false; // False on first step
      error = solve(isContinuation);
      if (error) {
         UTIL_THROW("Failure to converge initial state of sweep");
      } else {
         accept(sNew);
      }

      // Loop over states on path
      bool finished = false;   // Are we finished with the loop?
      while (!finished) {
         error = 1;
         while (error) {

            sNew = s(0) + ds; 
            std::cout << std::endl;
            std::cout << "Attempt s = " << sNew << std::endl;

            // Set a guess for adjustable state variables by continuation
            setGuess(sNew);

            // Set non-adjustable system parameters to new values
            setParameters(sNew);

            // Attempt iterative SCFT solution
            isContinuation = true;
            error = solve(isContinuation);

            if (error) {

               // Upon failure, reset state to last converged solution
               reset();

               // Decrease ds by half
               ds *= 0.50;
               if (ds < 0.1*ds0) {
                  UTIL_THROW("Step size too small in sweep");
               }

            } else {

               // Upon successful convergence, update history
               accept(sNew);

            }
         }
         if (sNew + ds > 1.0000001) {
            finished = true;
         }
      }

      // Clean up after end of sweep
      cleanup();

   }

   /*
   * Initialize history variables (must be called by setup function).
   */
   template <class State>
   void SweepTmpl<State>::initialize()
   {
      UTIL_CHECK(historyCapacity_  > 1);
      UTIL_CHECK(states_.capacity() == historyCapacity_);
      UTIL_CHECK(sHistory_.capacity() == historyCapacity_);
      UTIL_CHECK(stateHistory_.capacity() == historyCapacity_);

      // Check allocation of all states, allocate if necessary
      for (int i = 0; i < historyCapacity_; ++i) {
         checkAllocation(states_[i]);
      }

      // Set pointers in stateHistory_ to refer to objects in states_
      nAccept_ = 0;
      historySize_ = 0;
      for (int i = 0; i < historyCapacity_; ++i) {
         sHistory_[i] = 0.0;
         stateHistory_[i] = &states_[i];
      }
   }

   template <class State>
   void SweepTmpl<State>::accept(double sNew)
   {
      // Shift elements of sHistory_
      for (int i = historyCapacity_ - 1; i > 0; --i) {
         sHistory_[i] = sHistory_[i-1];
      }
      sHistory_[0] = sNew;

      // Shift elements of stateHistory_ (pointers to stored solutions)
      State* temp;
      temp = stateHistory_[historyCapacity_-1];
      for (int i = historyCapacity_ - 1; i > 0; --i) {
         stateHistory_[i] = stateHistory_[i-1];
      }
      stateHistory_[0] = temp;

      // Update counters
      ++nAccept_;
      if (historySize_ < historyCapacity_) {
         ++historySize_;
      }

      // Call getSolution to copy system state to state(0).
      getSolution();

   }

   template <class State>
   void SweepTmpl<State>::setCoefficients(double sNew)
   {
      if (historySize_ == 1) {
         c_[0] = 1.0;
      } else 
      if (historySize_ == 2) {
         c_[1] = -(sNew - s(0))/(sNew - s(1));
         c_[0] = 1.0 - c_[1];
      } else {
         UTIL_THROW("Invalid value of historySize");
      }
   }

   template <class State>
   void SweepTmpl<State>::cleanup()
   {}

} // namespace Pscf
#endif
