#ifndef RPG_LR_POST_AM_COMPRESSOR_TPP
#define RPG_LR_POST_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpg/System.h>
#include <rpg/fts/compressor/intra/IntraCorrelation.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg{

   using namespace Util;

   // Constructor
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
   : Compressor<D>(system),
     isAllocated_(false),
     intraCorrelation_(system)
   { setClassName("LrAmCompressor"); }

   // Destructor
   template <int D>
   LrAmCompressor<D>::~LrAmCompressor()
   {}

   // Read parameters from file
   template <int D>
   void LrAmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::readErrorType(in);
   }
   
   // Initialize just before entry to iterative loop.
   template <int D>
   void LrAmCompressor<D>::setup(bool isContinuation)
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::setup(isContinuation);
      
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
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         intraCorrelationK_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(meshSize);
            wFieldTmp_[i].allocate(meshSize);
         }
            
         isAllocated_ = true;
      }

      // Store value of initial guess chemical potential fields
      DArray<RField<D>> const * currSys = &system().w().rgrid();
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks,nThreads>>>(w0_[i].cField(), 
                                          (*currSys)[i].cField(), meshSize);
      }
      
      // Compute homopolymer intraCorrelation
      intraCorrelationK_ = intraCorrelation_.computeIntraCorrelations();
   }
   
   template <int D>
   int LrAmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::solve();
      //mdeCounter_ = AmIteratorTmpl<Compressor<D>,DArray<double>>::totalItr();
      return solve;
   }

   // Assign one array to another
   template <int D>
   void LrAmCompressor<D>::setEqual(Field<cudaReal>& a, Field<cudaReal> const & b)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cField(), b.cField(), a.capacity());
   }
   
   // Compute and return inner product of two vectors.
   template <int D>
   double LrAmCompressor<D>::dotProduct(Field<cudaReal> const & a, 
                                            Field<cudaReal> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = (double)gpuInnerProduct(a.cField(), b.cField(), n);
      return product;
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double LrAmCompressor<D>::maxAbs(Field<cudaReal> const & a)
   {
      int n = a.capacity();
      cudaReal max = gpuMaxAbs(a.cField(), n);
      return (double)max;
   }

   // Update basis
   template <int D>
   void 
   LrAmCompressor<D>::updateBasis(RingBuffer< Field<cudaReal> > & basis,
                                      RingBuffer< Field<cudaReal> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // New basis vector is difference between two most recent states
      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            (hists[0].cField(), hists[1].cField(), newBasis_.cField(),n);

      basis.append(newBasis_);
   }

   template <int D>
   void
   LrAmCompressor<D>::addHistories(Field<cudaReal>& trial,
                                       RingBuffer<Field<cudaReal> > const & basis,
                                       DArray<double> coeffs,
                                       int nHist)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(trial.capacity(), nBlocks, nThreads);

      for (int i = 0; i < nHist; i++) {
         pointWiseAddScale<<<nBlocks, nThreads>>>
            (trial.cField(), basis[i].cField(), -1*coeffs[i], trial.capacity());
      }
   }
   
   template<int D>
   double LrAmCompressor<D>::computeLambda(double r)
   {
      return 1.0;
   }

   template <int D>
   void LrAmCompressor<D>::addPredictedError(Field<cudaReal>& fieldTrial,
                                                 Field<cudaReal> const & resTrial,
                                                 double lambda)
   {
      int n = fieldTrial.capacity();
      const double vMonomer = system().mixture().vMonomer();
      const int meshSize = system().domain().mesh().size();
      
       // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);
      
      // Convert resTrial to RField<D> type
      assignReal<<<nBlocks, nThreads>>>(resid_.cField(), resTrial.cField(), meshSize);
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      
      // Combine with Linear response factor to update second step
      scaleComplex<<<nBlocks, nThreads>>>(residK_.cField(), 1.0/vMonomer, kSize_);
      inPlacePointwiseDivComplex<<<nBlocks, nThreads>>>(residK_.cField(), intraCorrelationK_.cField(), kSize_);
      
      // Convert back to real Space
      system().fft().inverseTransform(residK_, resid_);
      
      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cField(), resid_.cField(), lambda, fieldTrial.capacity());
   }

   // Does the system have an initial field guess?
   template <int D>
   bool LrAmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return the number of elements in a field vector
   template <int D>
   int LrAmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   // Get the current field from the system
   template <int D>
   void LrAmCompressor<D>::getCurrent(Field<cudaReal>& curr)
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
            ((*currSys)[0].cField(), w0_[0].cField(), curr.cField(), meshSize); 
   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void LrAmCompressor<D>::evaluate()
   {  
      system().compute(); 
      ++mdeCounter_;
   }

   // Compute the residual for the current system state
   template <int D>
   void LrAmCompressor<D>::getResidual(Field<cudaReal>& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      // Initialize residuals
      assignUniformReal<<<nBlocks, nThreads>>>(resid.cField(), -1.0, meshSize);

      // Compute SCF residual vector elements
      for (int i = 0; i < nMonomer; i++) {
         pointWiseAdd<<<nBlocks, nThreads>>>
            (resid.cField(), system().c().rgrid(i).cField(), meshSize);
      }

   }

   // Update the current system field coordinates
   template <int D>
   void LrAmCompressor<D>::update(Field<cudaReal>& newGuess)
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
            (w0_[i].cField(), newGuess.cField(), wFieldTmp_[i].cField(), meshSize);
      }
      
      // Set system r grid
      system().setWRGrid(wFieldTmp_);
   }
   
   template<int D>
   void LrAmCompressor<D>::outputToLog()
   {}
   
   template<int D>
   void LrAmCompressor<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "LrAmCompressor time contributions:\n";
      AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::outputTimers(out);
   }
   
   // Clear timers and MDE counter 
   template<int D>
   void LrAmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::clearTimers();
      mdeCounter_ = 0;
   }
   
}
}
#endif
