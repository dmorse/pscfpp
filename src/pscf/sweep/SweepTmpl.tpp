#ifndef PSCF_SWEEP_TMPL_TPP
#define PSCF_SWEEP_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepTmpl.h"
#include <util/misc/Log.h>

namespace Pscf {

   using namespace Util;

   /*
   * Constructor
   */
   template <class State>
   SweepTmpl<State>::SweepTmpl(int historyCapacity)
    : ns_(0),
      baseFileName_(),
      writeCRGrid_(false),
      writeCBasis_(false),
      writeWRGrid_(false),
      historyCapacity_(historyCapacity),
      historySize_(0),
      nAccept_(0),
      reuseState_(true)
   {  setClassName("SweepTmpl"); }

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
      // Default values
      baseFileName_ = "";

      // Read parameters
      read<int>(in, "ns", ns_);
      readOptional<std::string>(in, "baseFileName", baseFileName_);
      readOptional<int>(in, "historyCapacity", historyCapacity_);
      readOptional<bool>(in, "reuseState", reuseState_);
      readOptional<bool>(in, "writeCRGrid", writeCRGrid_);
      readOptional<bool>(in, "writeCBasis", writeCBasis_);
      readOptional<bool>(in, "writeWRGrid", writeWRGrid_);

      // Allocate required arrays
      UTIL_CHECK(historyCapacity_ > 0);
      states_.allocate(historyCapacity_);
      stateHistory_.allocate(historyCapacity_);
      sHistory_.allocate(historyCapacity_);
      c_.allocate(historyCapacity_);
   }

   template <class State>
   void SweepTmpl<State>::sweep()
   {

      // Compute and output ds
      double ds = 1.0/double(ns_);
      double ds0 = ds;
      Log::file() << std::endl;
      Log::file() << "ns = " << ns_ << std::endl;
      Log::file() << "ds = " << ds  << std::endl;

      // Initial setup, before a sweep
      setup();

      // Solve for initial state of sweep
      double sNew = 0.0;
      Log::file() << std::endl;
      Log::file() << "===========================================\n";
      Log::file() << "Attempt s = " << sNew << std::endl;

      int error;
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

         // Loop over iteration attempts, with decreasing ds as needed
         error = 1;
         while (error) {

            // Set a new contour variable value sNew
            sNew = s(0) + ds; 
            Log::file() << std::endl;
            Log::file() << "===========================================\n";
            Log::file() << "Attempt s = " << sNew << std::endl;

            // Set non-adjustable system parameters to new values
            setParameters(sNew);

            // Guess new state variables by polynomial extrapolation.
            // This function must both compute the extrapolation and
            // set initial guess values in the parent system.
            extrapolate(sNew);

            // Attempt iterative SCFT solution
            isContinuation = reuseState_;
            error = solve(isContinuation);

            // Process success or failure
            if (error) {
               Log::file() << "Backtrack and halve sweep step size:" 
                           << std::endl;

               // Upon failure, reset state to last converged solution
               reset();

               // Decrease ds by half
               ds *= 0.50;
               if (ds < 0.1*ds0) {
                  UTIL_THROW("Sweep decreased ds too many times.");
               }

            } else {

               // Upon successful convergence, update history and nAccept
               accept(sNew);

            }
         }
         if (sNew + ds > 1.0000001) {
            finished = true;
         }
      }
      Log::file() << "===========================================\n";

      // Clean up after end of the entire sweep
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

      // Set pointers in stateHistory_ to refer to objects in array states_
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

      // Update counters nAccept_ and historySize_
      ++nAccept_;
      if (historySize_ < historyCapacity_) {
         ++historySize_;
      }

      // Call getSolution to copy system state to state(0).
      getSolution();
   }

   /*
   * Use Lagrange polynomials to compute coefficients for continuation.
   */
   template <class State>
   void SweepTmpl<State>::setCoefficients(double sNew)
   {
      UTIL_CHECK(historySize_ <= historyCapacity_);
      if (historySize_ == 1) {
         c_[0] = 1.0;
      } else {
         double num, den;
         int i, j;
         for (i = 0; i < historySize_; ++i) {
            num = 1.0;
            den = 1.0;
            for (j = 0; j < historySize_; ++j) {
               if (j != i) {
                  num *= (sNew - s(j));
                  den *= (s(i) - s(j));
               }
            }
            c_[i] = num/den;
         }
      }
      // Theory: The coefficient c_[i] is a Lagrange polynomial 
      // function c_[i](sNew) of sNew that is defined such that
      // c_[i](s(i)) = 1 and c_[i](s(j)) = 0 for any j != i.
      // Given a set of values y(i) previously obtained at contour
      // variable values s(i) for all 0 <= i < historySize_, a 
      // linear combination f(sNew) = \sum_i c_[i](sNew)*y(i) summed 
      // over i = 0, ..., historySize_ - 1 gives a polynomial in
      // sNew that passes through all previous values, such that
      // f(sNew) = y(i) for sNew = s(i).
   }

   /*
   * Clean up after the end of a sweep (empty default implementation).
   */
   template <class State>
   void SweepTmpl<State>::cleanup()
   {}

} // namespace Pscf
#endif
