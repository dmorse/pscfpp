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
         newBasis_.allocate(meshSize);
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
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      // Store current fields
      DArray<RField<D>> const * currSys = &system().w().rgrid();
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks,nThreads>>>(w0_[i].cArray(), 
                                          (*currSys)[i].cArray(), meshSize);
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
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), a.capacity());
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

      const int n = hists[0].capacity();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // New basis vector is difference between two most recent states
      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            (hists[0].cArray(), hists[1].cArray(), newBasis_.cArray(),n);

      basis.append(newBasis_);

   }

   template <int D>
   void
   LrAmPreCompressor<D>::addHistories(DeviceArray<cudaReal>& trial,
                                   RingBuffer<DeviceArray<cudaReal> > const & basis,
                                   DArray<double> coeffs,
                                   int nHist)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(trial.capacity(), nBlocks, nThreads);

      for (int i = 0; i < nHist; i++) {
         pointWiseAddScale<<<nBlocks, nThreads>>>
            (trial.cArray(), basis[i].cArray(), -1*coeffs[i], trial.capacity());
      }
   }

   template <int D>
   void 
   LrAmPreCompressor<D>::addPredictedError(DeviceArray<cudaReal>& fieldTrial,
                                           DeviceArray<cudaReal> const & resTrial,
                                           double lambda)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);

      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cArray(), resTrial.cArray(), lambda, fieldTrial.capacity());
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
      // Straighten out fields into  linear arrays
      const int meshSize = system().domain().mesh().size();
      
      // Pointer to fields on system
      DArray<RField<D>> const * currSys = &system().w().rgrid();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      /*
      * The field that we are adjusting is the Langrange multiplier field 
      * with number of grid pts components.The current value is the difference 
      * between w and w0_ for the first monomer (any monomer should give the same answer)
      */
      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            ((*currSys)[0].cArray(), w0_[0].cArray(), curr.cArray(), meshSize); 
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
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vMonomer = system().mixture().vMonomer();  
      
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
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      
      // Residual combine with Linear response factor
      scaleComplex<<<nBlocks, nThreads>>>(residK_.cArray(), 1.0/vMonomer, kSize_);
      inPlacePointwiseDivComplex<<<nBlocks, nThreads>>>(residK_.cArray(), 
                                                        intraCorrelation_.cArray(), 
                                                        kSize_);
   
      // Convert back to real Space
      system().fft().inverseTransform(residK_, resid_);
      
      // Assign resid_ to resid
      assignReal<<<nBlocks, nThreads>>>(resid.cArray(), resid_.cArray(), meshSize);
      
   }

   // Update the current system field coordinates
   template <int D>
   void LrAmPreCompressor<D>::update(DeviceArray<cudaReal>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // New field is the w0_ + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++){
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>
            (w0_[i].cArray(), newGuess.cArray(), wFieldTmp_[i].cArray(), meshSize);
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
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Initialize residuals to -1
      assignUniformReal<<<nBlocks, nThreads>>>(error_.cArray(), -1.0, meshSize);
     
      // Add composition of each monomer
      for (int i = 0; i < nMonomer; i++) {
         pointWiseAdd<<<nBlocks, nThreads>>>
            (error_.cArray(), system().c().rgrid(i).cArray(), meshSize);
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
