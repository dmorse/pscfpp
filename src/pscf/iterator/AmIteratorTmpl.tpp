#ifndef PSCF_AM_ITERATOR_TMPL_TPP
#define PSCF_AM_ITERATOR_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NanException.h"
#include <pscf/inter/Interaction.h>
#include <pscf/math/LuSolver.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/misc/Timer.h>
#include <util/misc/FileMaster.h>

#include <cmath>

namespace Pscf
{

   using namespace Util;

   // Public member functions

   /*
   * Constructor
   */
   template <typename Iterator, typename T>
   AmIteratorTmpl<Iterator,T>::AmIteratorTmpl()
    : epsilon_(0),
      maxItr_(200),
      maxHist_(50),
      verbose_(1),
      errorType_("relNormResid"),
      error_(0.0),
      nBasis_(0),
      itr_(0),
      totalItr_(0),
      nElem_(0),
      lambda_(1.0),
      r_(0.9),
      useLambdaRamp_(true),
      isAllocatedAM_(false)
   {  ParamComposite::setClassName("AmIteratorTmpl"); }

   /*
   * Destructor
   */
   template <typename Iterator, typename T>
   AmIteratorTmpl<Iterator,T>::~AmIteratorTmpl()
   {}

   /*
   * Read parameter file block.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::readParameters(std::istream& in)
   {

      // Error tolerance for scalar error. Stop when error < epsilon.
      read(in, "epsilon", epsilon_);

      // Maximum number of iterations (optional)
      // Default set in constructor (200) 
      readOptional(in, "maxItr", maxItr_);

      // Maximum number of previous vectors in history (optional)
      // Default set in constructor (50) 
      readOptional(in, "maxHist", maxHist_);

      // Verbosity level of error reporting, values 0-3 (optional)
      // Initialized to 1 by default in constructor
      // verbose_ = 0 => FTS
      // verbose_ = 1 => concise
      // verbose_ = 2 => report all error measures at end
      // verbose_ = 3 => report all error measures every iteration
      readOptional(in, "verbose", verbose_);

      #if 0
      // Ramping parameter for correction step
      // Default set to 0.9. lambda = 1 - r_^nBasis for nBasis < maxHist
      readOptional(in, "r", r_);
      #endif

      #ifdef PSCF_AM_TEST
      // Default set to be false
      readOptional(in, "hasAmTest", hasAmTest_);
      #endif
   }

   /*
   * Iteratively solve a fixed-point problem.
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::solve(bool isContinuation)
   {
      // Initialization and allocate operations on entry to loop.
      setup(isContinuation);

      // Preconditions for generic Anderson-mixing (AM) algorithm.
      UTIL_CHECK(hasInitialGuess());
      UTIL_CHECK(isAllocatedAM_);

      // Start overall timer
      timerTotal_.start();

      // Solve MDE for initial state
      timerMDE_.start();
      evaluate();
      timerMDE_.stop();

      // Iterative loop
      nBasis_ = stateBasis_.size();
      for (itr_ = 0; itr_ < maxItr_; ++itr_) {

         // Append current state to stateHistory_ ringbuffer
         getCurrent(temp_);
         stateHistory_.append(temp_);

         timerAM_.start();

         if (verbose_ > 2) {
            Log::file() << "------------------------------- \n";
         }

         if (verbose_ > 0){
            Log::file() << " Iteration " << Int(itr_,5);
         }

         // Compute current residual vector, store in temp_
         timerResid_.start();
         getResidual(temp_);

         // Append current residual to residualHistory_ ring buffer
         residualHistory_.append(temp_);
         timerResid_.stop();

         // Compute scalar error used to test convergence
         timerError_.start();
         try {
            error_ = computeError(verbose_);
         } catch (const NanException&) {
            Log::file() << ",  error  =             NaN" << std::endl;
            break; // Exit loop if a NanException is caught
         }
         if (verbose_ > 0 && verbose_ < 3) {
            Log::file() << ",  error  = " << Dbl(error_, 15)
                        << std::endl;
         }
         timerError_.stop();

         // Output additional details of this iteration to the log file
         outputToLog();

         #ifdef PSCF_AM_TEST
         // Compute errors and ratios used for algorithm testing
         if (hasAmTest_){
            // Compute errors and ratios used for algorithm testing
            correctionError_ = computeError(0);
            if (itr_ > 0) {
               projectionRatio_ += projectionError_/preError_;
               correctionRatio_ += correctionError_/preError_;
               testCounter ++;
            }
            preError_ = correctionError_;
         }
         #endif

         // Check for convergence
         if (error_ < epsilon_) {

            // Stop timers
            timerAM_.stop();
            timerTotal_.stop();

            if (verbose_ > 2) {
               Log::file() << "-------------------------------\n";
            }

            if (verbose_ > 0) {
               Log::file() << " Converged\n";
            }

            // Output error report if not done previously
            if (verbose_ == 2) {
               Log::file() << "\n";
               computeError(2);
            }

            totalItr_ += itr_;

            // Successful completion (i.e., converged within tolerance)
            return 0;

         } else {

            // Compute optimal coefficients for basis vectors
            timerCoeff_.start();
            computeTrialCoeff();
            timerCoeff_.stop();

            // Compute updated trial state and residual vectors
            timerOmega_.start();
            updateTrial();

            // Correction (or "mixing") step
            addCorrection(stateTrial_, residualTrial_);

            // Update the parent system using new trial state
            update(stateTrial_);

            timerOmega_.stop();
            timerAM_.stop();

            // Perform the main calculation of the parent system -
            // Solve MDEs, compute phi's, compute stress if needed
            timerMDE_.start();
            evaluate();
            timerMDE_.stop();
         }

      }

      // Failure: Counter itr_ reached maxItr_ without converging
      timerTotal_.stop();
      Log::file() << "Iterator failed to converge.\n";
      return 1;

   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::outputTimers(std::ostream& out) const
   {
      // Output timing results, if requested.
      double total = timerTotal_.time();
      out << "\n";
      out << "                          ";
      out << "Total" << std::setw(22)<< "Per Iteration"
          << std::setw(9) << "Fraction" << "\n";
      out << "MDE solution:             "
          << Dbl(timerMDE_.time(), 9, 3)  << " s,  "
          << Dbl(timerMDE_.time()/totalItr_, 9, 3)  << " s,  "
          << Dbl(timerMDE_.time()/total, 9, 3) << "\n";
      out << "residual computation:     "
          << Dbl(timerResid_.time(), 9, 3)  << " s,  "
          << Dbl(timerResid_.time()/totalItr_, 9, 3)  << " s,  "
          << Dbl(timerResid_.time()/total, 9, 3) << "\n";
      out << "mixing coefficients:      "
          << Dbl(timerCoeff_.time(), 9, 3)  << " s,  "
          << Dbl(timerCoeff_.time()/totalItr_, 9, 3)  << " s,  "
          << Dbl(timerCoeff_.time()/total, 9, 3) << "\n";
      out << "checking convergence:     "
          << Dbl(timerError_.time(), 9, 3)  << " s,  "
          << Dbl(timerError_.time()/totalItr_, 9, 3)  << " s,  "
          << Dbl(timerError_.time()/total, 9, 3) << "\n";
      out << "updating guess:           "
          << Dbl(timerOmega_.time(), 9, 3)  << " s,  "
          << Dbl(timerOmega_.time()/totalItr_, 9, 3)  << " s,  "
          << Dbl(timerOmega_.time()/total, 9, 3)<< "\n";
      out << "total time:               "
          << Dbl(total, 9, 3) <<  " s,  "
          << Dbl(total/totalItr_, 9, 3) << " s  \n";
      out << "\n";

      #ifdef PSCF_AM_TEST
      if (hasAmTest_){
         out << "Average Projection Step Reduction Ratio:     "
             << Dbl(projectionRatio_/testCounter, 3, 3)<< "\n";
         out << "Average Correction Step Reduction Ratio:     "
             << Dbl(correctionRatio_/testCounter, 3, 3)<< "\n";
      }
      #endif
   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::clearTimers()
   {
      timerMDE_.clear();
      timerResid_.clear();
      timerError_.clear();
      timerCoeff_.clear();
      timerOmega_.clear();
      timerAM_.clear();
      timerTotal_.clear();
      totalItr_ = 0;
   }

   // Protected member functions - Initialization

   /*
   * Read and validate optional errorType string parameter.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::readErrorType(std::istream& in)
   {
      // Note: errorType_ is initialized to "relNormResid" in constructor

      // Read optional errorType_ string if present
      readOptional(in, "errorType", errorType_);

      // Check validity, normalize to standard form
      bool isValid = isValidErrorType();

      // If not valid, throw Exception with error message
      if (!isValid) {
         std::string msg = "Invalid iterator error type [";
         msg += errorType_;
         msg += "] in parameter file";
         UTIL_THROW(msg.c_str());
      }
   }

   /*
   * Check validity of errorType_ string and normalize if valid.
   */
   template <typename Iterator, typename T>
   bool AmIteratorTmpl<Iterator,T>::isValidErrorType()
   {
      // Process possible synonyms
      if (errorType_ == "norm") errorType_ = "normResid";
      if (errorType_ == "rms") errorType_ = "rmsResid";
      if (errorType_ == "max") errorType_ = "maxResid";
      if (errorType_ == "relNorm") errorType_ = "relNormResid";

      // Check value
      bool valid;
      valid = (errorType_ == "normResid"
            || errorType_ == "rmsResid"
            || errorType_ == "maxResid"
            || errorType_ == "relNormResid");
      return valid;

   }

   /*
   * Read and validate optional errorType string parameter.
   */
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::readMixingParameters(std::istream& in,
                                                    bool useLambdaRamp)
   {
      // Set default parameter for useLambdaRamp_
      useLambdaRamp_ = useLambdaRamp;

      // Optionally read mixing parameters
      readOptional(in, "lambda", lambda_); 
      readOptional(in, "useLambdaRamp", useLambdaRamp_); 
      if (useLambdaRamp_) {
         readOptional(in, "r", r_);
      }
   }

   /*
   * Allocate memory required by the AM algorithm, if needed.
   * (protected, non-virtual)
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::allocateAM()
   {
      // If already allocated, do nothing and return
      if (isAllocatedAM_) return;

      // Compute and set number of elements in a residual vector
      nElem_ = nElements();

      // Allocate ring buffers
      stateHistory_.allocate(maxHist_+1);
      residualHistory_.allocate(maxHist_+1);
      stateBasis_.allocate(maxHist_);
      residualBasis_.allocate(maxHist_);

      // Allocate arrays used in iteration
      stateTrial_.allocate(nElem_);
      residualTrial_.allocate(nElem_);
      temp_.allocate(nElem_);

      // Allocate arrays/matrices used in coefficient calculation
      U_.allocate(maxHist_, maxHist_);
      v_.allocate(maxHist_);
      coeffs_.allocate(maxHist_);

      isAllocatedAM_ = true;
   }

   /*
   * Initialize just before entry to iterative loop.
   * (virtual)
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setup(bool isContinuation)
   {
      if (!isAllocatedAM()) {
         allocateAM();
      } else {
         // Clear residual and state history buffers
         residualHistory_.clear();
         stateHistory_.clear();
         if (!isContinuation) {
            // Clear bases iff not a continuation
            residualBasis_.clear();
            stateBasis_.clear();
         }
      }
   }

   /*
   * Clear all history and basis vector data.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::clear()
   {
      UTIL_CHECK(isAllocatedAM_);
      residualHistory_.clear();
      stateHistory_.clear();
      residualBasis_.clear();
      stateBasis_.clear();
   }

   // Private non-virtual member functions

   /*
   * Compute optimal coefficients of basis vectors.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::computeTrialCoeff()
   {
      // If first iteration and bases are empty
      // then initialize U, v and coeff arrays
      if (itr_ == 0 && nBasis_ == 0) {
         int m, n;
         for (m = 0; m < maxHist_; ++m) {
            v_[m] = 0.0;
            coeffs_[m] = 0.0;
            for (n = 0; n < maxHist_; ++n) {
               U_(m, n) = 0.0;
            }
         }
         return;
      }

      // Do nothing else on first iteration
      if (itr_ == 0) return;

      // Update stateBasis_, residualBasis_, U_ matrix and v_ vector
      if (stateHistory_.size() > 1) {

         // Update basis spanning differences of past state vectors
         updateBasis(stateBasis_, stateHistory_);

         // Update basis spanning differences of past residual vectors
         updateBasis(residualBasis_, residualHistory_);

         // Update the U matrix and v vector.
         updateU(U_, residualBasis_);

      }

      // Update nBasis_
      nBasis_ = residualBasis_.size();
      UTIL_CHECK(nBasis_ > 0);
      UTIL_CHECK(stateBasis_.size() == nBasis_);

      // Update v_ vector (dot product of basis vectors and residual)
      computeV(v_, residualHistory_[0], residualBasis_);

      // Solution of the matrix problem U coeffs_ = -v yields a coeff
      // vector that minimizes the norm of the residual. Below, we:
      //    1. Solve matrix equation: U coeffs_ = v
      //    2. Flip the sign of the resulting coeffs_ vector

      if (nBasis_ == 1) {
         // Solve explicitly for coefficient
         coeffs_[0] = v_[0] / U_(0,0);
      } else
      if (nBasis_ < maxHist_) {

         // Create temporary smaller version of U_, v_, coeffs_ .
         // This avoids reallocation of U_ with each iteration.
         DMatrix<double> tempU;
         DArray<double> tempv, tempcoeffs;
         tempU.allocate(nBasis_,nBasis_);
         tempv.allocate(nBasis_);
         tempcoeffs.allocate(nBasis_);
         for (int i = 0; i < nBasis_; ++i) {
            tempv[i] = v_[i];
            for (int j = 0; j < nBasis_; ++j) {
               tempU(i,j) = U_(i,j);
            }
         }

         // Solve matrix equation
         LuSolver solver;
         solver.allocate(nBasis_);
         solver.computeLU(tempU);
         solver.solve(tempv, tempcoeffs);

         // Transfer solution to full-sized member variable
         for (int i = 0; i < nBasis_; ++i) {
            coeffs_[i] = tempcoeffs[i];
         }

      } else
      if (nBasis_ == maxHist_) {

         LuSolver solver;
         solver.allocate(maxHist_);
         solver.computeLU(U_);
         solver.solve(v_, coeffs_);

      }

      // Flip the sign of the coeffs_ vector (see comment above)
      for (int i = 0; i < nBasis_; ++i) {
         coeffs_[i] *= -1.0;
      }

      return;
   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::updateTrial()
   {

      // Set state and residual trial vectors to current values
      setEqual(stateTrial_, stateHistory_[0]);
      setEqual(residualTrial_, residualHistory_[0]);

      // Add linear combinations of state and residual basis vectors
      if (nBasis_ > 0) {

         // Combine basis vectors into trial guess and predicted residual
         addEqVectors(stateTrial_, stateBasis_, coeffs_);
         addEqVectors(residualTrial_, residualBasis_, coeffs_);

      }

      #ifdef PSCF_AM_TEST
      // Additional direct computation of error after projection step
      if (hasAmTest_){
         timerOmega_.stop();

         // Additional MDE solution
         timerMDE_.start();
         update(stateTrial_);
         evaluate();
         timerMDE_.stop();

         // Residual computation
         timerResid_.start();
         getResidual(temp_);
         timerResid_.stop();

         // Compute error after projection
         timerError_.start();
         projectionError_ = computeError(temp_, stateTrial_, errorType_, 0);
         timerError_.stop();
         timerOmega_.start();
      }
      #endif

   }

   // Private virtual member functions

   /*
   * Compute ramped mixing parameter for Anderson mixing algorithm.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::computeLambda()
   {
      if (useLambdaRamp_ && (nBasis_ < maxHist_)) {
         double factor = 1.0 - pow(r_, nBasis_ + 1);
         return factor*lambda_;
      } else {
         return lambda_;
      }
   }

   /*
   * Update entire U matrix.
   */
   template <typename Iterator, typename T>
   void
   AmIteratorTmpl<Iterator, T> ::updateU(
                                 DMatrix<double> & U,
                                 RingBuffer<T> const & residualBasis)
   {
      int nBasis = residualBasis.size();
      int maxHist = U.capacity1();
      UTIL_CHECK(maxHist >= nBasis);

      // Update matrix U by shifting elements diagonally
      for (int m = maxHist - 1; m > 0; --m) {
         for (int n = maxHist-1; n > 0; --n) {
            U(m,n) = U(m-1,n-1);
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nBasis; ++m) {
         double dotprod = dotProduct(residualBasis[0],residualBasis[m]);
         if (m == 0) {
            U(0,0) = dotprod;
         } else {
            U(m,0) = dotprod;
            U(0,m) = dotprod;
         }
      }

   }

   /* 
   * Compute v vector (dot products of residual and basis vectors).
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator, T>::computeV(
                                   DArray<double> & v,
                                   T const & resCurrent,
                                   RingBuffer<T> const & residualBasis)
   {
      int nBasis = residualBasis.size();
      UTIL_CHECK(v.capacity() >= nBasis);
      for (int m = 0; m < nBasis; ++m) {
         v[m] = dotProduct(resCurrent, residualBasis[m]);
      }
   }

   /*
   * Compute scalar error, possibly output to log.
   */
   template <typename Iterator, typename T>
   double 
   AmIteratorTmpl<Iterator,T>::computeError(T&residTrial, 
                                            T&stateTrial,
                                            std::string errorType,
                                            int verbose)
   {
      double error = 0.0;

      // Find max residual vector element
      double maxRes  = maxAbs(residTrial);

      // Find norm of residual vector
      double normRes = norm(residTrial);

      // Check if calculation has diverged (normRes will be NaN)
      UTIL_CHECK(!std::isnan(normRes));

      // Find root-mean-squared residual element value
      double rmsRes = normRes/sqrt(nElements());

      // Find norm of residual vector relative to state
      double normField = norm(stateTrial);
      double relNormRes = normRes/normField;

      // Set error value
      if (errorType == "maxResid") {
         error = maxRes;
      } else if (errorType == "normResid") {
         error = normRes;
      } else if (errorType == "rmsResid") {
         error = rmsRes;
      } else if (errorType == "relNormResid") {
         error = relNormRes;
      } else {
         UTIL_THROW("Invalid iterator error type in parameter file.");
      }

      if (verbose > 1) {
         Log::file() << "\n";
         Log::file() << "Max Residual  = " << Dbl(maxRes,15) << "\n";
         Log::file() << "Residual Norm = " << Dbl(normRes,15) << "\n";
         Log::file() << "RMS Residual  = " << Dbl(rmsRes,15) << "\n";
         Log::file() << "Relative Norm = " << Dbl(relNormRes,15)
                     << std::endl;
      }

      return error;
   }

   /*
   * Compute error (use current residual, errorType).
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::computeError(int verbose)
   {
      return computeError(residualHistory_[0], stateHistory_[0], 
                          errorType_, verbose);
   }


   /*
   * Update a list of basis vectors.
   */
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::updateBasis(RingBuffer<T> & basis, 
                                           RingBuffer<T> const & history)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(history.size() >= 2);

      // Set up array to store new basis
      basis.advance();
      if (basis[0].isAllocated()) {
         UTIL_CHECK(basis[0].capacity() == history[0].capacity());
      } else {
         basis[0].allocate(history[0].capacity());
      }

      // New basis vector is difference between two most recent vectors
      subVV(basis[0], history[0], history[1]);
   }

   /*
   * Add a linear combination of basis vectors to a vector v. 
   */
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::addEqVectors(T& v, 
                                            RingBuffer<T> const & basis, 
                                            DArray<double> coeffs)
   {
      int nBasis = basis.size();
      UTIL_CHECK(coeffs.capacity() >= nBasis);
      for (int i = 0; i < nBasis; i++) {
         addEqVc(v, basis[i], coeffs[i]);
      }
   }

   /*
   * Add a correction based on the predicted remaining residual. 
   */
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::addCorrection(T& stateTrial, 
                                             T const & residualTrial)
   {
      double lambda = computeLambda();
      addEqVc(stateTrial, residualTrial, lambda);
   }

   // Private virtual vector math functions

   /*
   * Compute L2 norm of a vector.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::norm(T const & a)
   {
      double normSq = dotProduct(a, a);
      return sqrt(normSq);
   }

}
#endif
