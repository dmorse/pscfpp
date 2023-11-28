#ifndef PSCF_AM_ITERATOR_TMPL_TPP
#define PSCF_AM_ITERATOR_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/inter/Interaction.h>
#include <pscf/math/LuSolver.h>
#include "NanException.h"
#include <util/containers/FArray.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/misc/Timer.h>
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
    : errorType_("relNormResid"),
      epsilon_(0),
      lambda_(0),
      maxItr_(200),
      maxHist_(50),
      nBasis_(0),
      totalItr_(0),
      nElem_(0),
      verbose_(1),
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
      // Default set in constructor (200) or before calling this function
      readOptional(in, "maxItr", maxItr_);

      // Maximum number of previous states in history (optional)
      // Default set in constructor (50) or before calling this function
      readOptional(in, "maxHist", maxHist_);

      // Verbosity level of error reporting, values 0-3 (optional)
      // Initialized to 1 by default in constructor
      // verbose_ = 0 => FTS
      // verbose_ = 1 => concise
      // verbose_ = 2 => report all error measures at end
      // verbose_ = 3 => report all error measures every iteration
      readOptional(in, "verbose", verbose_);
   }

   /*
   * Iteratively solve a SCF problem.
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


         if (nBasis_ < maxHist_) {
            lambda_ = 1.0 - pow(0.9, nBasis_ + 1);
         } else {
            lambda_ = 1.0;
         }

         // Compute residual vector
         timerResid_.start();
         getResidual(temp_);

         // Append current residual to resHists_ ringbuffer
         resHists_.append(temp_);
         timerResid_.stop();

         // Compute scalar error, output report to log file.
         timerError_.start();
         double error;
         try {
            error = computeError(verbose_);
         } catch (const NanException&) {
            Log::file() << ",  error  =             NaN" << std::endl;
            break; // Exit loop if a NanException is caught
         }
         if (verbose_ > 0 && verbose_ < 3) {
             Log::file() << ",  error  = " << Dbl(error, 15) << std::endl;
         }
         timerError_.stop();

         // Output additional details of this iteration to the log file
         outputToLog();

         // Check for convergence
         if (error < epsilon_) {

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

            // Compute the trial updated field and update the system
            timerOmega_.start();
            updateGuess();
            timerOmega_.stop();

            timerAM_.stop();

            // Perform the main calculation of the parent system -
            // solve MDEs, compute phi's, compute stress if needed
            timerMDE_.start();
            evaluate();
            timerMDE_.stop();

         }

      }
      
      
      // Failure: iteration counter itr reached maxItr without converging
      timerTotal_.stop();

      Log::file() << "Iterator failed to converge.\n";
      return 1;

   }

   // Protected member functions

   /*
   * Set value of maxItr.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setMaxItr(int maxItr)
   {  maxItr_ = maxItr; }

   /*
   * Set value of maxHist (number of retained previous states)
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setMaxHist(int maxHist)
   {  maxHist_ = maxHist; }

   /*
   * Set and validate value of error type string.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setErrorType(std::string errorType)
   {  
      errorType_ = errorType; 

      if (!isValidErrorType()) {
         std::string msg = "Invalid iterator error type [";
         msg += errorType;
         msg += "] input in AmIteratorTmpl::setErrorType";
         UTIL_THROW(msg.c_str());
      }
   }

   /*
   * Read and validate errorType string parameter.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::readErrorType(std::istream& in)
   {
      // Read optional string
      // Note: errorType_ is initialized to "relNormResid" in constructor
      readOptional(in, "errorType", errorType_);

      if (!isValidErrorType()) {
         std::string msg = "Invalid iterator error type [";
         msg += errorType_;
         msg += "] in parameter file";
         UTIL_THROW(msg.c_str());
      }

   }

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
   * Allocate memory required by the Anderson-Mixing algorithm, if needed.
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

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::clear()
   {
      if (!isAllocatedAM_) return;

      // Clear histories and bases (ring buffers)
      if (verbose_ > 0) {
         Log::file() << "Clearing Anderson-mixing history.\n";
      }
      resHists_.clear();
      fieldHists_.clear();
      resBasis_.clear();
      fieldBasis_.clear();

      return;
   }

   // Private non-virtual member functions

   /*
   * Compute coefficients of basis vectors.
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

      // Update basis spanning differences of past field vectors
      updateBasis(fieldBasis_, fieldHists_);

      // Update basis spanning differences of past residual vectors
      updateBasis(resBasis_, resHists_);

      // Update nBasis_
      nBasis_ = fieldBasis_.size();
      UTIL_CHECK(fieldBasis_.size() == nBasis_);

      // Update the U matrix and v vector.
      updateU(U_, resBasis_, nBasis_);
      updateV(v_, resHists_[0], resBasis_, nBasis_);
      // Note: resHists_[0] is the current residual vector

      // Solve matrix equation problem to compute coefficients
      // that minmize the L2 norm of the residual vector.
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
         solver.solve(tempv,tempcoeffs);

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

      // Correct for predicted error
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
         if (!isContinuation) {
            clear();
         }
      }
   }

   /*
   * Compute L2 norm of a vector.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::norm(T const & a)
   {
      double normSq = dotProduct(a, a);
      return sqrt(normSq);
   }

   // Update entire U matrix
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
         U(m,0) = dotprod;
         U(0,m) = dotprod;
      }
   }

   template <typename Iterator, typename T>
   void 
   AmIteratorTmpl<Iterator, T>::updateV(DArray<double> & v,
                                        T const & resCurrent,
                                        RingBuffer<T> const & resBasis,
                                        int nHist)
   {
      for (int m = 0; m < nHist; ++m) {
         v[m] = dotProduct(resCurrent, resBasis[m]);
      }
   }

   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::computeError(int verbose)
   {
      double error = 0.0;

      if (verbose > 1) {

         Log::file() << "\n";

         // Find max residual vector element
         double maxRes  = maxAbs(resHists_[0]);
         Log::file() << "Max Residual  = " << Dbl(maxRes,15) << "\n";
   
         // Find norm of residual vector
         double normRes = norm(resHists_[0]);
         Log::file() << "Residual Norm = " << Dbl(normRes,15) << "\n";
   
         // Find root-mean-squared residual element value
         double rmsRes = normRes/sqrt(nElem_);
         Log::file() << "RMS Residual  = " << Dbl(rmsRes,15) << "\n";

         // Find norm of residual vector relative to field
         double normField = norm(fieldHists_[0]);
         double relNormRes = normRes/normField;
         Log::file() << "Relative Norm = " << Dbl(relNormRes,15) << std::endl;
   
         // Check if calculation has diverged (normRes will be NaN)
         UTIL_CHECK(!std::isnan(normRes));
   
         // Set error value
         if (errorType_ == "maxResid") {
            error = maxRes;
         } else if (errorType_ == "normResid") {
            error = normRes;
         } else if (errorType_ == "rmsResid") {
            error = rmsRes;
         } else if (errorType_ == "relNormResid") {
            error = relNormRes;
         } else {
            UTIL_THROW("Invalid iterator error type in parameter file.");
         }

      } else {

         // Set error value
         if (errorType_ == "maxResid") {
            error = maxAbs(resHists_[0]);
         } else if (errorType_ == "normResid") {
            error = norm(resHists_[0]);
         } else if (errorType_ == "rmsResid") {
            error = norm(resHists_[0])/sqrt(nElem_);
         } else if (errorType_ == "relNormResid") {
            double normRes = norm(resHists_[0]);
            double normField = norm(fieldHists_[0]);
            error = normRes/normField;
         } else {
            UTIL_THROW("Invalid iterator error type in parameter file.");
         }
         //Log::file() << ",  error  = " << Dbl(error, 15) << "\n";

      }

      return error;
   }
   
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::outputTimers(std::ostream& out)
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

}
#endif
