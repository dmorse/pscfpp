#ifndef RPG_LR_POST_AM_COMPRESSOR_TPP
#define RPG_LR_POST_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpg/System.h>
#include <rpg/fts/compressor/intra/IntraCorrelation.h>
#include <prdc/cuda/resources.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg{

   using namespace Util;

   // Constructor
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
   : Compressor<D>(system),
     intra_(system),
     isIntraCalculated_(false),
     isAllocated_(false)
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
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::readErrorType(in);
   }
   
   // Initialize just before entry to iterative loop.
   template <int D>
   void LrAmCompressor<D>::setup(bool isContinuation)
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::setup(isContinuation);
      
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
      if (!isAllocated_) {
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         intraCorrelationK_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(dimensions);
            wFieldTmp_[i].allocate(dimensions);
         }
            
         isAllocated_ = true;
      }
      
      // Compute intraCorrelation
      if (!isIntraCalculated_){
         intra_.computeIntraCorrelations(intraCorrelationK_);
         isIntraCalculated_ = true;
      }

      // Store value of initial guess chemical potential fields
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0_[i], system().w().rgrid(i));
      }
      
   }
   
   template <int D>
   int LrAmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::solve();
      //mdeCounter_ = AmIteratorTmpl<Compressor<D>,DArray<double>>::totalItr();
      return solve;
   }

   // Assign one array to another
   template <int D>
   void LrAmCompressor<D>::setEqual(DeviceArray<cudaReal>& a, 
                                    DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b); 
   }
   
   // Compute and return inner product of two vectors.
   template <int D>
   double LrAmCompressor<D>::dotProduct(DeviceArray<cudaReal> const & a, 
                                        DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b); 
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double LrAmCompressor<D>::maxAbs(DeviceArray<cudaReal> const & a)
   {
      return Reduce::maxAbs(a);
   }

   // Update basis
   template <int D>
   void 
   LrAmCompressor<D>::updateBasis(RingBuffer< DeviceArray<cudaReal> > & basis,
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
   LrAmCompressor<D>::addHistories(DeviceArray<cudaReal>& trial,
                                   RingBuffer<DeviceArray<cudaReal> > const & basis,
                                   DArray<double> coeffs,
                                   int nHist)
   {
      for (int i = 0; i < nHist; i++) {
         VecOp::addEqVc(trial, basis[i], -1.0 * coeffs[i]);
      }
   }
   
   template<int D>
   double LrAmCompressor<D>::computeLambda(double r)
   {
      return 1.0;
   }

   template <int D>
   void LrAmCompressor<D>::addPredictedError(DeviceArray<cudaReal>& fieldTrial,
                                             DeviceArray<cudaReal> const & resTrial,
                                             double lambda)
   {
      int n = fieldTrial.capacity();
      const double vMonomer = system().mixture().vMonomer();
      
      // Convert resTrial to RField<D> type
      VecOp::eqV(resid_, resTrial);
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      
      // Combine with Linear response factor to update second step
      VecOp::divEqVc(residK_, intraCorrelationK_, vMonomer);
      
      // Convert back to real space (destroys residK_)
      system().fft().inverseTransformUnsafe(residK_, resid_);
      
      // fieldTrial += resid_ * lambda
      VecOp::addEqVc(fieldTrial, resid_, lambda);
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
   void LrAmCompressor<D>::getCurrent(DeviceArray<cudaReal>& curr)
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
   void LrAmCompressor<D>::evaluate()
   {  
      system().compute(); 
      ++mdeCounter_;
   }

   // Compute the residual for the current system state
   template <int D>
   void LrAmCompressor<D>::getResidual(DeviceArray<cudaReal>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();

      // Initialize resid to c field of species 0 minus 1
      VecOp::subVS(resid, system().c().rgrid(0), 1.0);

      // Add other c fields to get SCF residual vector elements
      for (int i = 1; i < nMonomer; i++) {
         VecOp::addEqV(resid, system().c().rgrid(i));
      }
   }

   // Update the current system field coordinates
   template <int D>
   void LrAmCompressor<D>::update(DeviceArray<cudaReal>& newGuess)
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
   void LrAmCompressor<D>::outputToLog()
   {}
   
   template<int D>
   void LrAmCompressor<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "LrAmCompressor time contributions:\n";
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::outputTimers(out);
   }
   
   // Clear timers and MDE counter 
   template<int D>
   void LrAmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >::clearTimers();
      mdeCounter_ = 0;
   }
   
}
}
#endif
