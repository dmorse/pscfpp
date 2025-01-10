#ifndef RPG_LR_AM_COMPRESSOR_TPP
#define RPG_LR_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmPreCompressor.h"
#include <rpg/System.h>
#include <rpg/fts/compressor/intra/IntraCorrelation.h>  
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/chem/Monomer.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/cuda/GpuResources.h>
#include <complex>
#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Constructor
   template <int D>
   LrAmPreCompressor<D>::LrAmPreCompressor(System<D>& system)
    : Compressor<D>(system),
      isAllocated_(false),
      intra_(system)
   {  setClassName("LrAmPreCompressor"); }

   // Destructor
   template <int D>
   LrAmPreCompressor<D>::~LrAmPreCompressor()
   {}

   // Read parameters from file
   template <int D>
   void LrAmPreCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::readErrorType(in);
   }
      
   // Initialize just before entry to iterative loop.
   template <int D>
   void LrAmPreCompressor<D>::setup(bool isContinuation)
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, 
                     DeviceArray<cudaReal> >::setup(isContinuation);
      
      // Compute Fourier space kMeshDimensions_
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
      
      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_){
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         error_.allocate(dimensions);
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         intraCorrelation_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(dimensions);
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }
      
      // Store value of initial guess chemical potential fields
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0_[i], system().w().rgrid(i));
      }
      
      // Compute intramolecular correlation
      intraCorrelation_ = intra_.computeIntraCorrelations();
   }
  
   // Iterative solver (AM algorithm) 
   template <int D>
   int LrAmPreCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, 
                                 DeviceArray<cudaReal> >::solve();
      //mdeCounter_ = AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal>>::totalItr();
      return solve;
   }

   // Assign one array to another
   template <int D>
   void LrAmPreCompressor<D>::setEqual(DeviceArray<cudaReal>& a, 
                                       DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b); 
   }

   // Compute and return inner product of two vectors.
   template <int D>
   double LrAmPreCompressor<D>::dotProduct(DeviceArray<cudaReal> const & a, 
                                           DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double LrAmPreCompressor<D>::maxAbs(DeviceArray<cudaReal> const & a)
   {
      return Reduce::maxAbs(a);
   }

   // Update basis
   template <int D>
   void 
   LrAmPreCompressor<D>::updateBasis(RingBuffer< DeviceArray<cudaReal> > & basis,
                               RingBuffer< DeviceArray<cudaReal> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      basis.advance();
      if (basis[0].isAllocated()) {
         UTIL_CHECK(basis[0].capacity() == hists[0].capacity());
      } else {
         basis[0].allocate(hists[0].capacity());
      }

      VecOp::subVV(basis[0], hists[0], hists[1]);
   }

   template <int D>
   void
   LrAmPreCompressor<D>::addHistories(DeviceArray<cudaReal>& trial,
                                   RingBuffer<DeviceArray<cudaReal> > const & basis,
                                   DArray<double> coeffs,
                                   int nHist)
   {
      for (int i = 0; i < nHist; i++) {
         VecOp::addEqVc(trial, basis[i], -1.0 * coeffs[i]);
      }
   }

   template <int D>
   void 
   LrAmPreCompressor<D>::addPredictedError(DeviceArray<cudaReal>& fieldTrial,
                                           DeviceArray<cudaReal> const & resTrial,
                                           double lambda)
   {
      VecOp::addEqVc(fieldTrial, resTrial, lambda);
   }

   // Does the system have an initial field guess?
   template <int D>
   bool LrAmPreCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return the number of elements in a field vector
   template <int D>
   int LrAmPreCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   // Get the current field from the system
   template <int D>
   void LrAmPreCompressor<D>::getCurrent(DeviceArray<cudaReal>& curr)
   {
      /*
      * The field that we are adjusting is the Langrange multiplier field 
      * with number of grid pts components.The current value is the difference 
      * between w and w0_ for the first monomer (any monomer should give the 
      * same answer)
      */
      VecOp::subVV(curr, system().w().rgrid(0), w0_[0]); 
   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void LrAmPreCompressor<D>::evaluate()
   {  
      system().compute(); 
      ++mdeCounter_;
   }

   // Compute the residual for the current system state
   template <int D>
   void LrAmPreCompressor<D>::getResidual(DeviceArray<cudaReal>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();  

      // Compute incompressibility constraint error vector elements
      VecOp::subVS(resid_, system().c().rgrid(0), 1.0);
      for (int i = 1; i < nMonomer; i++) {
         VecOp::addEqV(resid_, system().c().rgrid(i));
      }
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      
      // Residual combine with Linear response factor
      VecOp::divEqVc(residK_, intraCorrelation_, vMonomer);
   
      // Convert back to real Space
      system().fft().inverseTransform(residK_, resid_);
      
      // Assign resid_ to resid
      VecOp::eqV(resid, resid_);
      
   }

   // Update the current system field coordinates
   template <int D>
   void LrAmPreCompressor<D>::update(DeviceArray<cudaReal>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      
      // New field is the w0_ + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(wFieldTmp_[i], w0_[i], newGuess);
      }
      
      // Set system r grid
      system().setWRGrid(wFieldTmp_);
   }

   template<int D>
   void LrAmPreCompressor<D>::outputToLog()
   {}
   
   template<int D>
   void LrAmPreCompressor<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Compressor times contributions:\n";
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::outputTimers(out);
   }
   
   
   template<int D>
   void LrAmPreCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::clearTimers();
      mdeCounter_ = 0;
   }
      
   template<int D>
   double LrAmPreCompressor<D>::computeLambda(double r)
   {
      return 1.0;
   }
   
   template<int D>
   double LrAmPreCompressor<D>::computeError(DeviceArray<cudaReal>& residTrial, 
                                          DeviceArray<cudaReal>& fieldTrial,
                                          std::string errorType,
                                          int verbose)
   {
      double error = 0.0;
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      
      // Initialize residuals to sum of all monomer compositions minus 1
      VecOp::subVS(error_, system().c().rgrid(0), 1.0);
      for (int i = 1; i < nMonomer; i++) {
         VecOp::addEqV(error_, system().c().rgrid(i));
      }
      
      // Find max residual vector element
      double maxRes  = maxAbs(error_);
     
      // Find norm of residual vector
      double normRes = AmIteratorTmpl<Compressor<D>, 
                                      DeviceArray<cudaReal> >::norm(error_);
      
      // Find root-mean-squared residual element value
      double rmsRes = normRes/sqrt(n);
      if (verbose > 1) {
         Log::file() << "\n";
         Log::file() << "Max Residual  = " << Dbl(maxRes,15) << "\n";
         Log::file() << "Residual Norm = " << Dbl(normRes,15) << "\n";
         Log::file() << "RMS Residual  = " << Dbl(rmsRes,15);

         // Check if calculation has diverged (normRes will be NaN)
         UTIL_CHECK(!std::isnan(normRes));
         error = normRes;
         Log::file() <<"\n";
         // Set error value
         if (errorType == "maxResid") {
            error = maxRes;
         } else if (errorType == "normResid") {
            error = normRes;
         } else if (errorType == "rmsResid") {
            error = rmsRes;
         } else {
            UTIL_THROW("Invalid iterator error type in parameter file.");
         }

      } else {
         // Set error value
         if (errorType == "maxResid") {
            error = maxRes;
         } else if (errorType == "normResid") {
            error = normRes;
         } else if (errorType == "rmsResid") {
            error = normRes/sqrt(n);
         } else {
            UTIL_THROW("Invalid iterator error type in parameter file.");
         }
         //Log::file() << ",  error  = " << Dbl(error, 15) << "\n";
      }

      return error;
   }
   
}
}
#endif
