#ifndef PSPG_AM_COMPRESSOR_TPP
#define PSPG_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <pspg/System.h>
#include <pspg/field/RDField.h>
#include <util/global.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;

   // Constructor
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
   : Compressor<D>(system),
     isAllocated_(false)
   { setClassName("AmCompressor"); }

   // Destructor
   template <int D>
   AmCompressor<D>::~AmCompressor()
   {  setClassName("AmCompressor"); }

   // Read parameters from file
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, DField<cudaReal> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, DField<cudaReal> >::readErrorType(in);
   }
   
      
   // Initialize just before entry to iterative loop.
   template <int D>
   void AmCompressor<D>::setup(bool isContinuation)
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DField<cudaReal> >::setup(isContinuation);
      
      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_){
         newBasis_.allocate(meshSize);
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(meshSize);
            wFieldTmp_[i].allocate(meshSize);
         }
         isAllocated_ = true;
      }
      
      // Store value of initial guess chemical potential fields
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      // Pointer to fields on system
      DArray<RDField<D>> const * currSys = &system().w().rgrid();
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks,nThreads>>>(w0_[i].cDField(), 
                                          (*currSys)[i].cDField(), meshSize);
         
      }
   }
   
   template <int D>
   int AmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DField<cudaReal> >::solve();
      //counter_ += AmIteratorTmpl<Compressor<D>, DField<cudaReal>>::totalItr(); 
      return solve;
   }

   // Assign one array to another
   template <int D>
   void AmCompressor<D>::setEqual(DField<cudaReal>& a, DField<cudaReal> const & b)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cDField(), b.cDField(), a.capacity());
   }

   // Compute and return inner product of two vectors.
   template <int D>
   double AmCompressor<D>::dotProduct(DField<cudaReal> const & a, 
                                      DField<cudaReal> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = (double)gpuInnerProduct(a.cDField(), b.cDField(), n);
      return product;
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double AmCompressor<D>::maxAbs(DField<cudaReal> const & a)
   {
      int n = a.capacity();
      cudaReal max = gpuMaxAbs(a.cDField(), n);
      return (double)max;
   }

   // Update basis
   template <int D>
   void 
   AmCompressor<D>::updateBasis(RingBuffer< DField<cudaReal> > & basis,
                                RingBuffer< DField<cudaReal> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      
      // GPU resources
      // New basis vector is difference between two most recent states
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);
      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            (hists[0].cDField(), hists[1].cDField(), newBasis_.cDField(),n);

      basis.append(newBasis_);

   }

   template <int D>
   void
   AmCompressor<D>::addHistories(DField<cudaReal>& trial,
                                 RingBuffer<DField<cudaReal> > const & basis,
                                 DArray<double> coeffs,
                                 int nHist)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(trial.capacity(), nBlocks, nThreads);

      for (int i = 0; i < nHist; i++) {
         pointWiseAddScale<<<nBlocks, nThreads>>>
            (trial.cDField(), basis[i].cDField(), -1*coeffs[i], trial.capacity());
      }
   }

   template <int D>
   void AmCompressor<D>::addPredictedError(DField<cudaReal>& fieldTrial,
                                           DField<cudaReal> const & resTrial,
                                           double lambda)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);

      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cDField(), resTrial.cDField(), lambda, fieldTrial.capacity());
   }

   // Does the system have an initial field guess?
   template <int D>
   bool AmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return the number of elements in a field vector
   template <int D>
   int AmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   // Get the current field from the system
   template <int D>
   void AmCompressor<D>::getCurrent(DField<cudaReal>& curr)
   {
      // Straighten out fields into  linear arrays
      const int meshSize = system().domain().mesh().size();
   
      // Pointer to fields on system
      DArray<RDField<D>> const * currSys = &system().w().rgrid();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      /*
      * The field that we are adjusting is the Langrange multiplier field 
      * with number of grid pts components.The current value is the difference 
      * between w and w0_ for the first monomer (any monomer should give the same answer)
      */
      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            ((*currSys)[0].cDField(), w0_[0].cDField(), curr.cDField(), meshSize);  

   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void AmCompressor<D>::evaluate()
   {  system().compute(); }

   // Compute the residual for the current system state
   template <int D>
   void AmCompressor<D>::getResidual(DField<cudaReal>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Initialize residuals to -1
      assignUniformReal<<<nBlocks, nThreads>>>(resid.cDField(), -1, meshSize);
      for (int i = 0; i < nMonomer; i++) {
         pointWiseAdd<<<nBlocks, nThreads>>>
            (resid.cDField(), system().c().rgrid(i).cDField(), meshSize);
      }
   }

   // Update the current system field coordinates
   template <int D>
   void AmCompressor<D>::update(DField<cudaReal>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      //New field is the w0_ + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++){
         pointWiseBinaryAdd<<<nBlocks, nThreads>>>
            (w0_[i].cDField(), newGuess.cDField(), wFieldTmp_[i].cDField(), meshSize);
      }
      
      // set system r grid
      system().setWRGrid(wFieldTmp_);
   }

   template<int D>
   void AmCompressor<D>::outputToLog()
   {}

}
}
#endif