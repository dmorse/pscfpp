#ifndef RPG_LR_COMPRESSOR_TPP
#define RPG_LR_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrCompressor.h"
#include <rpg/System.h>
#include <rpg/fts/compressor/intra/IntraCorrelation.h> 
#include <util/global.h>
#include <util/format/Dbl.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/iterator/NanException.h>
#include <prdc/crystal/shiftToMinimum.h>


namespace Pscf {
namespace Rpg{

   using namespace Util;

   // Constructor
   template <int D>
   LrCompressor<D>::LrCompressor(System<D>& system)
    : Compressor<D>(system),
      epsilon_(0.0),
      itr_(0),
      maxItr_(0),
      totalItr_(0),
      errorType_("rmsResid"),
      verbose_(0),
      isAllocated_(false),
      intraCorrelation_(system)
   { setClassName("LrCompressor"); }

   // Destructor
   template <int D>
   LrCompressor<D>::~LrCompressor()
   {}

   // Read parameters from file
   template <int D>
   void LrCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      maxItr_ = 60;
      read(in, "epsilon", epsilon_);
      readOptional(in, "maxItr", maxItr_);
      readOptional(in, "verbose", verbose_);
      readOptional(in, "errorType", errorType_);
   }

   // Initialize just before entry to iterative loop.
   template <int D>
   void LrCompressor<D>::setup()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
         }
      }
      
      // Compute number of points in k-space grid
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         kSize_ *= kMeshDimensions_[i];
      }
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Allocate memory required by AM algorithm if not done earlier.
      if (!isAllocated_){
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         wFieldTmp_.allocate(nMonomer);
         intraCorrelationK_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }
      
      // Compute intraCorrelation (homopolymer)
      intraCorrelationK_ = intraCorrelation_.computeIntraCorrelations();
   }

   template <int D>
   int LrCompressor<D>::compress()
   {
      // Initialization and allocate operations on entry to loop.
      setup();
      UTIL_CHECK(isAllocated_);

      // Start overall timer
      timerTotal_.start();

      // Solve MDE
      timerMDE_.start();
      system().compute();
      ++mdeCounter_;
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
         getResidual();
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
            ++mdeCounter_;
            timerMDE_.stop();

         }

      }

      // Failure: iteration counter itr reached maxItr without converging
      timerTotal_.stop();

      Log::file() << "Iterator failed to converge.\n";
      return 1;

   }

   // Compute the residual for the current system state
   template <int D>
   void LrCompressor<D>::getResidual()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      // Initialize residuals to -1
      assignUniformReal<<<nBlocks, nThreads>>>(resid_.cArray(), -1, meshSize);

      // Compute incompressibility constraint error vector elements
      for (int i = 0; i < nMonomer; i++) {
         pointWiseAdd<<<nBlocks, nThreads>>>
            (resid_.cArray(), system().c().rgrid(i).cArray(), meshSize);
      }
      
   }

   // update system w field using linear response approximation
   template <int D>
   void LrCompressor<D>::updateWFields(){
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vMonomer = system().mixture().vMonomer();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      
      // Compute change in fields using estimated Jacobian
      scaleComplex<<<nBlocks, nThreads>>>(residK_.cArray(), 1.0/vMonomer, kSize_);
      inPlacePointwiseDivComplex<<<nBlocks, nThreads>>>(residK_.cArray(), 
                                                        intraCorrelationK_.cArray(), 
                                                        kSize_);
   
      // Convert back to real Space
      system().fft().inverseTransform(residK_, resid_);
      
      // Update new fields
      for (int i = 0; i < nMonomer; i++){
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>
            (system().w().rgrid(i).cArray(), resid_.cArray(), 
             wFieldTmp_[i].cArray(), meshSize);
      }
      
      // Set system r grid
      system().setWRGrid(wFieldTmp_);
   }

   template<int D>
   void LrCompressor<D>::outputToLog()
   {}

   template<int D>
   void LrCompressor<D>::outputTimers(std::ostream& out)
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

   template<int D>
   void LrCompressor<D>::clearTimers()
   {
      timerTotal_.clear();
      timerMDE_.clear();
      mdeCounter_ = 0;
      totalItr_ = 0;
   }

   template <int D>
   double LrCompressor<D>::maxAbs(RField<D> const & a)
   {
      int n = a.capacity();
      cudaReal max = gpuMaxAbs(a.cArray(), n);
      return (double)max;
   }

   // Compute and return inner product of two vectors.
   template <int D>
   double LrCompressor<D>::dotProduct(RField<D> const & a,
                                      RField<D> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = (double)gpuInnerProduct(a.cArray(), b.cArray(), n);
      return product;
   }

   // Compute L2 norm of an RField
   template <int D>
   double LrCompressor<D>::norm(RField<D> const & a)
   {
      double normSq = dotProduct(a, a);
      return sqrt(normSq);
   }

   // Compute and return the scalar error
   template <int D>
   double LrCompressor<D>::computeError(int verbose)
   {
      const int meshSize = system().domain().mesh().size();
      double error = 0;

      // Find max residual vector element
      double maxRes  = maxAbs(resid_);
      // Find norm of residual vector
      double normRes = norm(resid_);
      // Find root-mean-squared residual element value
      double rmsRes = normRes/sqrt(meshSize);
      if (errorType_ == "maxResid") {
         error = maxRes;
      } else if (errorType_ == "normResid") {
         error = normRes;
      } else if (errorType_ == "rmsResid") {
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
