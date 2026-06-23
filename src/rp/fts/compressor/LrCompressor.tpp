#ifndef RP_LR_COMPRESSOR_TPP
#define RP_LR_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrCompressor.h"
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/iterator/NanException.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   LrCompressor<D,T>::LrCompressor(System<D,T>& system)
    : CompressorT(system),
      intra_(system),
      errorType_("rms"),
      epsilon_(0.0),
      itr_(0),
      maxItr_(100),
      totalItr_(0),
      verbose_(0),
      isAllocated_(false),
      isIntraCalculated_(false)
   {  ParamComposite::setClassName("LrCompressor"); }

   /*
   * Read parameters from file.
   */
   template <int D, class T>
   void LrCompressor<D,T>::readParameters(std::istream& in)
   {
      // Default values
      maxItr_ = 100;
      verbose_ = 0;
      errorType_ = "rms";

      // Read parameters
      ParamComposite::read(in, "epsilon", epsilon_);
      ParamComposite::readOptional(in, "maxItr", maxItr_);
      ParamComposite::readOptional(in, "verbose", verbose_);
      ParamComposite::readOptional(in, "errorType", errorType_);
   }

   /*
   * Initialize just before entry to iterative loop.
   */
   template <int D, class T>
   void LrCompressor<D,T>::setup()
   {
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      int kSize;
      FFTT::computeKMesh(dimensions, kMeshDimensions_, kSize);

      // Allocate memory required by AM algorithm if not done earlier.
      if (!isAllocated_) {
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         const int nMonomer = system().mixture().nMonomer();
         wFieldTmp_.allocate(nMonomer);
         intraCorrelationK_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }

      // Compute intraCorrelation
      if (!isIntraCalculated_){
         intra_.computeOmegaTotal(intraCorrelationK_);
         isIntraCalculated_ = true;
      }
   }

   /*
   * Adjust pressure field to find partial saddle point.
   */
   template <int D, class T>
   int LrCompressor<D,T>::compress()
   {
      // Initialization and allocate operations on entry to loop.
      setup();
      UTIL_CHECK(isAllocated_);

      // Start overall timer
      timerTotal_.start();

      // Solve MDE
      timerMDE_.start();
      system().compute();
      ++CompressorT::mdeCounter_;
      timerMDE_.stop();

      // Iterative loop
      for (itr_ = 0; itr_ < maxItr_; ++itr_) {

         if (verbose_ > 2) {
            Log::file() << "------------------------------- \n";
         }

         if (verbose_ > 0){
            Log::file() << " Iteration " << Int(itr_,5);
         }

         // Compute residual vector
         computeResidual();
         double error;
         try {
            error = computeError(verbose_);
         } catch (const NanException&) {
            Log::file() << ",  error  =             NaN" << std::endl;
            break; // Exit loop if a NanException is caught
         }
         if (verbose_ > 0) {
            Log::file() << ",  error  = " << Dbl(error, 15) << std::endl;
         }

         // Check for convergence
         if (error < epsilon_) {

            // Successful completion (i.e., converged within tolerance)
            timerTotal_.stop();

            if (verbose_ > 2) {
               Log::file() << "-------------------------------\n";
            }
            if (verbose_ > 0) {
               Log::file() << " Converged\n";
            }
            if (verbose_ == 2) {
               Log::file() << "\n";
               computeError(2);
            }
            //mdeCounter_ += itr_;
            totalItr_ += itr_;

            return 0; // Success

         } else{

            // Not yet converged.
            updateWFields();
            timerMDE_.start();
            system().compute();
            ++CompressorT::mdeCounter_;
            timerMDE_.stop();

         }

      }

      // Failure: iteration counter itr reached maxItr without converging
      timerTotal_.stop();
      Log::file() << "Iterator failed to converge.\n";
      return 1;

   }

   /*
   * Compute and store the residual for the current system state.
   */
   template <int D, class T>
   void LrCompressor<D,T>::computeResidual()
   {
      // Initialize resid_ to -1
      VecOp::eqS(resid_, -1.0);

      // Add all c fields
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::addEqV(resid_, system().c().rgrid(i));
      }
   }

   /*
   * Update system w fields using linear response approximation.
   */
   template <int D, class T>
   void LrCompressor<D,T>::updateWFields()
   {
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();

      // Convert residual to Fourier Space
      system().domain().fft().forwardTransform(resid_, residK_);

      // Compute change in fields using estimated Jacobian
      VecOp::divEqVc(residK_, intraCorrelationK_, vMonomer);

      // Convert back to real space (destroys residK_)
      system().domain().fft().inverseTransformUnsafe(residK_, resid_);

      // Update new fields
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(wFieldTmp_[i], system().w().rgrid(i), resid_);
      }

      // Set system w fields
      system().w().setRGrid(wFieldTmp_);
   }

   /*
   * Output information to a log file (do-nothing implementation).
   */
   template <int D, class T>
   void LrCompressor<D,T>::outputToLog()
   {}

   /*
   * Output timing information to evaluate performance.
   */
   template <int D, class T>
   void LrCompressor<D,T>::outputTimers(std::ostream& out) const
   {
      // Output timing results, if requested.
      double total = timerTotal_.time();
      out << "\n";
      out << "LrCompressor time contributions:\n";
      out << "                          ";
      out << "Total" << std::setw(22)<< "Per Iteration"
          << std::setw(9) << "Fraction" << "\n";
      out << "MDE solution:             "
          << Dbl(timerMDE_.time(), 9, 3)  << " s,  "
          << Dbl(timerMDE_.time()/totalItr_, 9, 3)  << " s,  "
          << Dbl(timerMDE_.time()/total, 9, 3) << "\n";
      out << "\n";
   }

   /*
   * Clear all timers.
   */
   template <int D, class T>
   void LrCompressor<D,T>::clearTimers()
   {
      timerTotal_.clear();
      timerMDE_.clear();
      CompressorT::mdeCounter_ = 0;
      totalItr_ = 0;
   }

   /*
   * Compute and return the scalar error.
   */
   template <int D, class T>
   double LrCompressor<D,T>::computeError(int verbose)
   {
      const int meshSize = system().domain().mesh().size();
      double error = 0.0;

      // Find max residual vector element
      double maxRes  = Reduce::maxAbs(resid_);
      // Find norm of residual vector
      double normRes = sqrt(Reduce::sumSq(resid_));
      // Find root-mean-squared residual element value
      double rmsRes = normRes/sqrt(meshSize);
      if (errorType_ == "max") {
         error = maxRes;
      } else if (errorType_ == "norm") {
         error = normRes;
      } else if (errorType_ == "rms") {
         error = rmsRes;
      } else {
         UTIL_THROW("Invalid iterator error type in parameter file.");
      }

      if (verbose > 1) {
         Log::file() << "\n";
         Log::file() << "Max Residual  = " << Dbl(maxRes,15) << "\n";
         Log::file() << "Residual Norm = " << Dbl(normRes,15) << "\n";
         Log::file() << "RMS Residual  = " << Dbl(rmsRes,15) << "\n";

         // Check if calculation has diverged (normRes will be NaN)
         UTIL_CHECK(!std::isnan(normRes));
      }
      return error;
   }

}
}
#endif
