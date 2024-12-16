#ifndef RPG_AM_COMPRESSOR_TPP
#define RPG_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpg/System.h>
#include <prdc/cuda/RField.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg{

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
   {}

   // Read parameters from file
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, 
                     DeviceArray<cudaReal> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, 
                     DeviceArray<cudaReal> >::readErrorType(in);
   }
   
      
   // Initialize just before entry to iterative loop.
   template <int D>
   void AmCompressor<D>::setup(bool isContinuation)
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const IntVec<D> dimensions = system().domain().mesh().dimensions();
      
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, 
                     DeviceArray<cudaReal> >::setup(isContinuation);
      
      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_){
         newBasis_.allocate(meshSize);
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
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
      // Pointer to fields on system
      DArray<RField<D>> const * currSys = &system().w().rgrid();
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks,nThreads>>>(w0_[i].cArray(), 
                                          (*currSys)[i].cArray(), meshSize);
         
      }
   }
   
   template <int D>
   int AmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, 
                                 DeviceArray<cudaReal> >::solve();
      //mdeCounter_ = AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal>>::totalItr(); 
      return solve;
   }

   // Assign one array to another
   template <int D>
   void AmCompressor<D>::setEqual(DeviceArray<cudaReal>& a, 
                                  DeviceArray<cudaReal> const & b)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), 
                                        a.capacity());
   }

   // Compute and return inner product of two vectors.
   template <int D>
   double AmCompressor<D>::dotProduct(DeviceArray<cudaReal> const & a, 
                                      DeviceArray<cudaReal> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = (double)gpuInnerProduct(a.cArray(), b.cArray(), n);
      return product;
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double AmCompressor<D>::maxAbs(DeviceArray<cudaReal> const & a)
   {
      int n = a.capacity();
      cudaReal max = gpuMaxAbs(a.cArray(), n);
      return (double)max;
   }

   // Update basis
   template <int D>
   void 
   AmCompressor<D>::updateBasis(RingBuffer< DeviceArray<cudaReal> > & basis,
                                RingBuffer< DeviceArray<cudaReal> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      
      // GPU resources
      // New basis vector is difference between two most recent states
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);
      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            (hists[0].cArray(), hists[1].cArray(), newBasis_.cArray(),n);

      basis.append(newBasis_);

   }

   template <int D>
   void
   AmCompressor<D>::addHistories(DeviceArray<cudaReal>& trial,
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
   void AmCompressor<D>::addPredictedError(DeviceArray<cudaReal>& fieldTrial,
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
   bool AmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return the number of elements in a field vector
   template <int D>
   int AmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   // Get the current field from the system
   template <int D>
   void AmCompressor<D>::getCurrent(DeviceArray<cudaReal>& curr)
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
   void AmCompressor<D>::evaluate()
   {  
      system().compute();
      ++mdeCounter_; 
   }

   // Compute the residual for the current system state
   template <int D>
   void AmCompressor<D>::getResidual(DeviceArray<cudaReal>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Initialize residuals to -1
      assignUniformReal<<<nBlocks, nThreads>>>(resid.cArray(), -1, meshSize);
      for (int i = 0; i < nMonomer; i++) {
         pointWiseAdd<<<nBlocks, nThreads>>>
            (resid.cArray(), system().c().rgrid(i).cArray(), meshSize);
      }
   }

   // Update the current system field coordinates
   template <int D>
   void AmCompressor<D>::update(DeviceArray<cudaReal>& newGuess)
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
            (w0_[i].cArray(), newGuess.cArray(), wFieldTmp_[i].cArray(), meshSize);
      }
      
      // set system r grid
      system().setWRGrid(wFieldTmp_);
   }
   
   template<int D>
   double AmCompressor<D>::computeLambda(double r)
   {
      return 1.0;
   }

   template<int D>
   void AmCompressor<D>::outputToLog()
   {}
   
   template<int D>
   void AmCompressor<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Compressor times contributions:\n";
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::outputTimers(out);
   }
   
   template<int D>
   void AmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::clearTimers();
      mdeCounter_ = 0;
   }

}
}
#endif