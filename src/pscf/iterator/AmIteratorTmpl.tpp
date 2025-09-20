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
      lambda_(0),
      r_(0.9),
      nBasis_(0),
      itr_(0),
      totalItr_(0),
      nElem_(0),
      isAllocatedAM_(false)
   {  setClassName("AmIteratorTmpl"); }

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

      // Maximum number of iterations (optional parameter)
      // Default set in constructor (200) 
      readOptional(in, "maxItr", maxItr_);

      // Maximum number of previous states in history (optional)
      // Default set in constructor (50) 
      readOptional(in, "maxHist", maxHist_);

      // Verbosity level of error reporting, values 0-3 (optional)
      // Initialized to 1 by default in constructor
      // verbose_ = 0 => FTS
      // verbose_ = 1 => concise
      // verbose_ = 2 => report all error measures at end
      // verbose_ = 3 => report all error measures every iteration
      readOptional(in, "verbose", verbose_);

      // Ramping parameter for correction step
      // Default set to 0.9. lambda = 1 - r_^Nh for Nh < maxHist
      readOptional(in, "correctionRamp", r_);

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
      nBasis_ = fieldBasis_.size();
      for (itr_ = 0; itr_ < maxItr_; ++itr_) {

         // Append current field to fieldHists_ ringbuffer
         getCurrent(temp_);
         fieldHists_.append(temp_);

         timerAM_.start();

         if (verbose_ > 2) {
            Log::file() << "------------------------------- \n";
         }

         if (verbose_ > 0){
            Log::file() << " Iteration " << Int(itr_,5);
         }

         lambda_ = computeLambda(r_);

         // Compute current residual vector, store in temp_
         timerResid_.start();
         getResidual(temp_);

         // Append current residual to resHists_ ring buffer
         resHists_.append(temp_);
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
            computeResidCoeff();
            timerCoeff_.stop();

            // Compute updated field and update the system
            timerOmega_.start();
            updateGuess();
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

   // Protected member functions

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
      fieldHists_.allocate(maxHist_+1);
      resHists_.allocate(maxHist_+1);
      fieldBasis_.allocate(maxHist_);
      resBasis_.allocate(maxHist_);

      // Allocate arrays used in iteration
      fieldTrial_.allocate(nElem_);
      resTrial_.allocate(nElem_);
      temp_.allocate(nElem_);

      // Allocate arrays/matrices used in coefficient calculation
      U_.allocate(maxHist_, maxHist_);
      v_.allocate(maxHist_);
      coeffs_.allocate(maxHist_);

      isAllocatedAM_ = true;
   }

   /*
   * Clear all history and basis vector data.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::clear()
   {
      if (!isAllocatedAM_) return;

      // Clear histories and bases (ring buffers)
      if (verbose_ > 0) {
         Log::file() << "Clearing AM field history and basis vectors.\n";
      }
      resHists_.clear();
      fieldHists_.clear();
      resBasis_.clear();
      fieldBasis_.clear();

      return;
   }

   // Private non-virtual member functions

   /*
   * Compute optimal coefficients of basis vectors.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::computeResidCoeff()
   {
      // If first iteration and history is empty
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
      }

      // Do nothing else on first iteration
      if (itr_ == 0) return;

      // Update fieldBasis_, resBasis_, U_ matrix and v_ vector
      if (fieldHists_.size() > 1) {

         // Update basis spanning differences of past field vectors
         updateBasis(fieldBasis_, fieldHists_);

         // Update basis spanning differences of past residual vectors
         updateBasis(resBasis_, resHists_);

         // Update nBasis_
         nBasis_ = fieldBasis_.size();
         UTIL_CHECK(fieldBasis_.size() == nBasis_);

         // Update the U matrix and v vector.
         updateU(U_, resBasis_, nBasis_);
      }

      UTIL_CHECK(nBasis_ > 0);

      // Update v_ vector (dot product of basis vectors and residual)
      computeV(v_, resHists_[0], resBasis_, nBasis_);

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
         // This is done to avoid reallocating U_ with each iteration.
         DMatrix<double> tempU;
         DArray<double> tempv,tempcoeffs;
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
   void AmIteratorTmpl<Iterator,T>::updateGuess()
   {

      // Set field and residual trial vectors to current values
      setEqual(fieldTrial_, fieldHists_[0]);
      setEqual(resTrial_, resHists_[0]);

      // Add linear combinations of field and residual basis vectors
      if (nBasis_ > 0) {

         // Combine basis vectors into trial guess and predicted residual
         addHistories(fieldTrial_, fieldBasis_, coeffs_, nBasis_);
         addHistories(resTrial_, resBasis_, coeffs_, nBasis_);

      }

      #ifdef PSCF_AM_TEST
      // Additional direct computation of error after projection step
      if (hasAmTest_){
         timerOmega_.stop();

         // Additional MDE solution
         timerMDE_.start();
         update(fieldTrial_);
         evaluate();
         timerMDE_.stop();

         // Residual computation
         timerResid_.start();
         getResidual(temp_);
         timerResid_.stop();

         // Compute error after projection
         timerError_.start();
         projectionError_ = computeError(temp_, fieldTrial_, errorType_, 0);
         timerError_.stop();
         timerOmega_.start();
      }
      #endif

      // Correction (or "mixing") step
      addPredictedError(fieldTrial_, resTrial_,lambda_);

      // Update system using new trial field
      update(fieldTrial_);


      return;
   }

   // Private virtual member functions

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
         // Clear residual and field history buffers
         resHists_.clear();
         fieldHists_.clear();
         if (!isContinuation) {
            // Clear bases iff not a continuation
            resBasis_.clear();
            fieldBasis_.clear();
         }
      }
   }

   /**
   * Compute ramped mixing parameter for Anderson mixing algorithm.
   * (virtual)
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::computeLambda(double r)
   {
      double lambda;
      if (nBasis_ < maxHist_) {
         lambda = 1.0 - pow(r, nBasis_ + 1);
      } else {
         lambda = 1.0;
      }
      return lambda;
   }

   /*
   * Update entire U matrix.
   */
   template <typename Iterator, typename T>
   void
   AmIteratorTmpl<Iterator, T> ::updateU(DMatrix<double> & U,
                                         RingBuffer<T> const & resBasis,
                                         int nHist)
   {
      // Update matrix U by shifting elements diagonally
      int maxHist = U.capacity1();
      for (int m = maxHist-1; m > 0; --m) {
         for (int n = maxHist-1; n > 0; --n) {
            U(m,n) = U(m-1,n-1);
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nHist; ++m) {
         double dotprod = dotProduct(resBasis[0],resBasis[m]);
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
   void
   AmIteratorTmpl<Iterator, T>::computeV(
                                   DArray<double> & v,
                                   T const & resCurrent,
                                   RingBuffer<T> const & resBasis,
                                   int nHist)
   {
      for (int m = 0; m < nHist; ++m) {
         v[m] = dotProduct(resCurrent, resBasis[m]);
      }
   }

   template <typename Iterator, typename T>
   double 
   AmIteratorTmpl<Iterator,T>::computeError(T&residTrial, 
                                            T&fieldTrial,
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

      // Find norm of residual vector relative to field
      double normField = norm(fieldTrial);
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
      return computeError(resHists_[0], fieldHists_[0], 
                          errorType_, verbose);
   }


   // Update the series of residual vectors.
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::updateBasis(RingBuffer<T> & basis, 
                                           RingBuffer<T> const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      // Set up array to store new basis
      basis.advance();
      if (basis[0].isAllocated()) {
         UTIL_CHECK(basis[0].capacity() == hists[0].capacity());
      } else {
         basis[0].allocate(hists[0].capacity());
      }

      // New basis vector is difference between two most recent states
      subVV(basis[0], hists[0], hists[1]);
   }

   /*
   * Compute trial field so as to minimize L2 norm of residual.
   */
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::addHistories(T& trial, 
                                            RingBuffer<T> const & basis, 
                                            DArray<double> coeffs, 
                                            int nHist)
   {
      for (int i = 0; i < nHist; i++) {
         addEqVc(trial, basis[i], coeffs[i]);
      }
   }

   /*
   * Add predicted error to the trial field.
   */
   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator,T>::addPredictedError(T& fieldTrial, 
                                                 T const & resTrial, 
                                                 double lambda)
   {
      addEqVc(fieldTrial, resTrial, lambda);
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

   #if 0
   /* 
   * Assign one vector to another: a = b.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setEqual(T& a, T const & b)
   {  a = b; }

   /*
   * Compute and return inner product of two vectors 
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::dotProduct(T const & a, 
                                                 T const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      double product = 0.0;
      for (int i = 0; i < n; i++) {
         // if either value is NaN, throw NanException
         if (std::isnan(a[i]) || std::isnan(b[i])) { 
            throw NanException("AmIteratorTmpl<Iterator,T>::dotProduct",
                               __FILE__,__LINE__,0);
         }
         product += a[i] * b[i];
      }
      return product;
   }

   /*
   * Compute and return maximum element of residual vector.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::maxAbs(T const & a)
   {
      const int n = a.capacity();
      double max = 0.0;
      double value;
      for (int i = 0; i < n; i++) {
         value = a[i];
         if (std::isnan(value)) { // if value is NaN, throw NanException
            throw NanException("AmIteratorTmpl<Iterator,T>::dotProduct",
                                __FILE__,__LINE__,0);
         }
         if (fabs(value) > max)
            max = fabs(value);
      }
      return max;
   }

   /*
   * Compute difference of two vectors.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::subVV(T& a,
                                          T const & b,
                                          T const & c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      UTIL_CHECK(n == c.capacity());
      for (int i = 0; i < n; i++) {
         a[i] = b[i] - c[i];
      }
   }

   /*
   * Compute and return inner product of two vectors 
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::addEqVc(T& a,
                                            T const & b,
                                            double c)
   {
      const int n = a.capacity();
      UTIL_CHECK(n == b.capacity());
      for (int i = 0; i < n; i++) {
         a[i] += c*b[i];
      }
   }
   #endif

}
#endif
