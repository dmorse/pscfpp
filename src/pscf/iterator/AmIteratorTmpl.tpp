#ifndef PSCF_AM_ITERATOR_TMPL_TPP
#define PSCF_AM_ITERATOR_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/inter/Interaction.h>
#include <util/containers/FArray.h>
#include <util/format/Dbl.h>
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
      nHist_(0),
      maxHist_(0),
      nElem_(0),
      diverged_(false),
      isAllocated_(false)
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
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);

      maxHist_ = 50;
      readOptional(in, "maxHist", maxHist_);
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

      if (!(errorType_ == "normResid"
         || errorType_ == "maxResid"
         || errorType_ == "relNormResid")) {
         UTIL_THROW("Invalid iterator error type in parameter file.");
      }
   }

   /*
   * Solve iteratively.
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::solve(bool isContinuation)
   {
      // Note: Parameter isContinuation is currently unused 

      // Initialization and allocate operations on entry to loop.
      setup();

      // Preconditions for generic algorithm.
      UTIL_CHECK(hasInitialGuess());
      UTIL_CHECK(isAllocated_);

      // Timers for analyzing performance
      Timer timerMDE;
      Timer timerStress;
      Timer timerAM;
      Timer timerResid;
      Timer timerConverged;
      Timer timerCoeff;
      Timer timerOmega;
      Timer timerTotal;

      // Start overall timer
      timerTotal.start();

      // Solve MDE for initial state
      timerMDE.start();
      evaluate();
      timerMDE.stop();

      // Iterative loop
      bool done;
      for (int itr = 0; itr < maxItr_; ++itr) {

         // Append current field to fieldHists_ ringbuffer
         getCurrent(temp_);
         fieldHists_.append(temp_);

         timerAM.start();

         Log::file()<<"---------------------"<<std::endl;
         Log::file()<<" Iteration  "<<itr<<std::endl;

         if (itr < maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr+1);
            nHist_ = itr;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         timerResid.start();
         // Compute residual vector
         getResidual(temp_);

         // Append current residual to resHists_ ringbuffer
         resHists_.append(temp_);
         timerResid.stop();

         // Test for convergence.
         // Also outputs error measures to Log::file()
         timerConverged.start();
         done = isConverged();
         timerConverged.stop();

         // Output additional details of this iteration to the log file
         outputToLog();

         if (diverged_) {

            // If calculation diverged, some residual values are
            // NaN. In this case, the calculation failed, and there
            // is no need to continue to the next iteration.
            break;
            
         } else if (done) {

            // Stop timers
            timerAM.stop();
            timerTotal.stop();

            Log::file() << "----------CONVERGED----------"<< std::endl;

            // Output timing results
            Log::file() << "\n";
            Log::file() << "Iterator times contributions:\n";
            Log::file() << "\n";
            Log::file() << "MDE solution:         "
                        << timerMDE.time()  << " s,  "
                        << timerMDE.time()/timerTotal.time() << "\n";
            Log::file() << "residual computation: "
                        << timerResid.time()  << " s,  "
                        << timerResid.time()/timerTotal.time() << "\n";
            Log::file() << "mixing coefficients:  "
                        << timerCoeff.time()  << " s,  "
                        << timerCoeff.time()/timerTotal.time() << "\n";
            Log::file() << "checking convergence: "
                        << timerConverged.time()  << " s,  "
                        << timerConverged.time()/timerTotal.time() << "\n";
            Log::file() << "updating guess:       "
                        << timerOmega.time()  << " s,  "
                        << timerOmega.time()/timerTotal.time() << "\n";
            Log::file() << "total time:           "
                        << timerTotal.time()   << " s  ";
            Log::file() << "\n\n";

            // Successful completion (i.e., converged within tolerance)
            cleanUp();
            return 0;

         } else {

            // Compute optimal coefficients for basis vectors
            timerCoeff.start();
            computeResidCoeff();
            timerCoeff.stop();

            // Compute the trial updated field and update the system
            timerOmega.start();
            updateGuess();
            timerOmega.stop();

            timerAM.stop();

            // Perform the main calculation of the parent system -
            // solve MDEs, compute phi's, compute stress if needed
            timerMDE.start();
            evaluate();
            timerMDE.stop();

         }

      }

      // Failure: iteration counter itr reached maxItr without converging
      timerTotal.stop();

      Log::file() << "Iterator failed to converge.\n";
      cleanUp();
      return 1;

   }

   // Protected member function

   /*
   * Allocate memory required by the Anderson-Mixing algorithm, if needed.
   * (protected, non-virtual)
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::allocateAM()
   {
      // If already allocated, do nothing and return
      if (isAllocated_) return;

      // Determine number of elements in a residual vector
      nElem_ = nElements();

      // Allocate ring buffers
      fieldHists_.allocate(maxHist_+1);
      fieldBasis_.allocate(maxHist_+1);
      resHists_.allocate(maxHist_+1);
      resBasis_.allocate(maxHist_+1);

      // Allocate arrays used in iteration
      fieldTrial_.allocate(nElem_);
      resTrial_.allocate(nElem_);
      temp_.allocate(nElem_);

      // Allocate arrays/matrices used in coefficient calculation
      U_.allocate(maxHist_, maxHist_);
      v_.allocate(maxHist_);
      coeffs_.allocate(maxHist_);

      isAllocated_ = true;
   }

   // Private non-virtual member functions

   // Private virtual member functions

   /*
   * Initialize just before entry to iterative loop.
   * (virtual)
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setup()
   {
      // Allocate memory required by AM algorithm if not done earlier.
      // If already allocated, do nothing and return.
      allocateAM();
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
   bool AmIteratorTmpl<Iterator,T>::isConverged()
   {

      // Find max residual vector element
      double maxRes  = maxAbs(resHists_[0]);
      Log::file() << "Max Residual  = " << Dbl(maxRes,15) << std::endl;

      // Find norm of residual vector
      double normRes = norm(resHists_[0]);
      Log::file() << "Residual Norm = " << Dbl(normRes,15) << std::endl;

      // Find norm of residual vector relative to field
      double relNormRes = normRes/norm(fieldHists_[0]);
      Log::file() << "Relative Norm = " << Dbl(relNormRes,15) << std::endl;

      // Check if calculation has diverged (normRes will be NaN)
      if (std::isnan(normRes)) {
         diverged_ = true;
         return false;
      }

      // Check if total error is below tolerance
      if (errorType_ == "normResid")
         return normRes < epsilon_;
      else if (errorType_ == "relNormResid")
         return relNormRes < epsilon_;
      else if (errorType_ == "maxResid")
         return maxRes < epsilon_;
      else
         UTIL_THROW("Invalid iterator error type in parameter file.");

   }

   /*
   * Compute coefficients of basis vectors.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::computeResidCoeff()
   {
      // Initialize matrix and vector of residual dot products
      // if this is the first iteration
      if (nHist_ == 0) {
         for (int m = 0; m < maxHist_; ++m) {
            v_[m] = 0.0;
            coeffs_[m] = 0.0;
            for (int n = 0; n < maxHist_; ++n) {
               U_(m, n) = 0.0;
            }
         }
         return;
      }

      // Update basis spanning differences of past residual vectors
      updateBasis(resBasis_, resHists_);

      // Update the U matrix and v vector.
      updateU(U_, resBasis_, nHist_);
      updateV(v_, resHists_[0], resBasis_, nHist_);
      // Note: resHists_[0] is the current residual vector

      // Solve matrix equation problem to compute coefficients
      // that minmize the L2 norm of the residual vector.
      if (nHist_ == 1) {
         // Solve explicitly for coefficient
         coeffs_[0] = v_[0] / U_(0,0);
      } else
      if (nHist_ < maxHist_) {

         // Create temporary smaller version of U_, v_, coeffs_ .
         // This is done to avoid reallocating U_ with each iteration.
         DMatrix<double> tempU;
         DArray<double> tempv,tempcoeffs;
         tempU.allocate(nHist_,nHist_);
         tempv.allocate(nHist_);
         tempcoeffs.allocate(nHist_);
         for (int i = 0; i < nHist_; ++i) {
            tempv[i] = v_[i];
            for (int j = 0; j < nHist_; ++j) {
               tempU(i,j) = U_(i,j);
            }
         }

         // Solve matrix equation
         LuSolver solver;
         solver.allocate(nHist_);
         solver.computeLU(tempU);
         solver.solve(tempv,tempcoeffs);

         // Transfer solution to full-sized member variable
         for (int i = 0; i < nHist_; ++i) {
            coeffs_[i] = tempcoeffs[i];
         }
      } else
      if (nHist_ == maxHist_) {
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
      if (nHist_ > 0) {

         // Update basis spanning differences of past field vectors
         updateBasis(fieldBasis_,fieldHists_);

         // Combine basis vectors into trial guess and predicted residual
         addHistories(fieldTrial_, fieldBasis_, coeffs_, nHist_);
         addHistories(resTrial_, resBasis_, coeffs_, nHist_);

      }

      // Correct for predicted error
      addPredictedError(fieldTrial_, resTrial_,lambda_);

      // Update system using new trial field
      update(fieldTrial_);

      return;
   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::cleanUp()
   {
      // Clear ring buffers
      resHists_.clear();
      resBasis_.clear();
      fieldHists_.clear();
      fieldBasis_.clear();
      return;
   }

}
#endif
