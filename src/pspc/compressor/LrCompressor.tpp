#ifndef PSPC_LR_COMPRESSOR_TPP
#define PSPC_LR_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrCompressor.h"
#include <pspc/System.h>
#include <util/global.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/iterator/NanException.h>
#include <prdc/crystal/shiftToMinimum.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   // Constructor
   template <int D>
   LrCompressor<D>::LrCompressor(System<D>& system)
   : Compressor<D>(system),
     counter_(0),
     epsilon_(0.0),
     itr_(0),
     maxItr_(0),
     errorType_("rmsResid"),
     verbose_(0),
     isAllocated_(false)
   {  setClassName("LrCompressor"); }

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
      // Allocate memory required by AM algorithm if not done earlier.
      if (!isAllocated_){
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         intraCorrelation_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(meshSize);
            wFieldTmp_[i].allocate(meshSize);
         }
         isAllocated_ = true;
      }
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j< meshSize; ++j){
            w0_[i][j] = system().w().rgrid(i)[j];
         }
      }
      computeIntraCorrelation();
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
      timerMDE_.stop();
      
      // Iterative loop
      for (itr_ = 0; itr_ < maxItr_; ++itr_) {
         if (verbose_ > 2) {
            Log::file() << "------------------------------- \n";
         }
         
         if (verbose_ > 0){
            Log::file() <<  std::endl;
            Log::file() << " Iteration " << Int(itr_,5) << std::endl;
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
         
         // Check for convergence
         if (error < epsilon_) {
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
            counter_ += itr_;
            // Successful completion (i.e., converged within tolerance)
            return 0;
         } else{
            updateWFields();
            timerMDE_.start();
            system().compute();
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

      // Initialize residuals
      for (int i = 0 ; i < meshSize; ++i) {
         resid_[i] = -1.0;
      }

       // Compute SCF residual vector elements
      for (int j = 0; j < nMonomer; ++j) {
        for (int k = 0; k < meshSize; ++k) {
           resid_[k] += system().c().rgrid(j)[k];
        }
      }
   }
   
   // update system w field using linear response approximation
   template <int D>
   void LrCompressor<D>::updateWFields(){  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size(); 
      const double vMonomer = system().mixture().vMonomer();   
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      MeshIterator<D> iter;
      iter.setDimensions(residK_.dftDimensions());
      for (iter.begin(); !iter.atEnd(); ++iter) {
         residK_[iter.rank()][0] *= 1.0 / (vMonomer * intraCorrelation_[iter.rank()]);
         residK_[iter.rank()][1] *= 1.0 / (vMonomer * intraCorrelation_[iter.rank()]);
      }
      
      // Convert back to real Space
      system().fft().inverseTransform(residK_, resid_);
      
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            wFieldTmp_[i][k] = system().w().rgrid(i)[k] + resid_[k];
         }
      }
      system().setWRGrid(wFieldTmp_);
   }

   template<int D>
   void LrCompressor<D>::outputToLog()
   {}
   
   template<int D>
   void LrCompressor<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Lr Compressor times contributions:\n";
   }
   
   template<int D>
   double LrCompressor<D>::computeDebye(double x)
   {
      if (x == 0){
         return 1.0;
      } else {
         return 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
   }
   
   template<int D>
   double LrCompressor<D>::computeIntraCorrelation(double qSquare)
   {
      const int np = system().mixture().nPolymer();
      const double vMonomer = system().mixture().vMonomer();
      // Overall intramolecular correlation
      double omega = 0;
      int monomerId; int nBlock; 
      double kuhn; double length; double g; double rg2; 
      Polymer<D> const * polymerPtr;
      
      for (int i = 0; i < np; i++){
         polymerPtr = &system().mixture().polymer(i);
         nBlock = polymerPtr->nBlock();
         double totalLength = 0;
         for (int j = 0; j < nBlock; j++) {
            totalLength += polymerPtr->block(j).length();
         }
         for (int j = 0; j < nBlock; j++) {
            monomerId = polymerPtr-> block(j).monomerId();
            kuhn = system().mixture().monomer(monomerId).kuhn();
            // Get the length (number of monomers) in this block.
            length = polymerPtr-> block(j).length();
            rg2 = length * kuhn* kuhn /6.0;
            g = computeDebye(qSquare*rg2);
            omega += length/totalLength * length * g/ vMonomer;
         }
      }
      return omega;
   }
   
   template<int D>
   void LrCompressor<D>::computeIntraCorrelation()
   {
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
         Gsq = system().unitCell().ksq(Gmin);
         intraCorrelation_[iter.rank()] = computeIntraCorrelation(Gsq);
      }
   }
   
      
   template<int D>
   void LrCompressor<D>::clearTimers()
   {

   }
   
   template <int D>
   double LrCompressor<D>::maxAbs(RField<D> const & a)
   {
      const int n = a.capacity();
      double max = 0.0;
      double value;
      for (int i = 0; i < n; i++) {
         value = a[i];
         if (std::isnan(value)) { // if value is NaN, throw NanException
            throw NanException("LrCompressor::dotProduct",__FILE__,__LINE__,0);
         }
         if (fabs(value) > max) {
            max = fabs(value);
         }
      }
      return max;
   }
   
   // Compute and return inner product of two vectors.
   template <int D>
   double LrCompressor<D>::dotProduct(RField<D> const & a, 
                                      RField<D> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = 0.0;
      for (int i = 0; i < n; i++) {
         // if either value is NaN, throw NanException
         if (std::isnan(a[i]) || std::isnan(b[i])) { 
            throw NanException("AmCompressor::dotProduct",__FILE__,__LINE__,0);
         }
         product += a[i] * b[i];
      }
      return product;
   }
   
   // Compute L2 norm 
   template <int D>
   double LrCompressor<D>::norm(RField<D> const & a)
   {
      double normSq = dotProduct(a, a);
      return sqrt(normSq);
   }

   
   template<int D>
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