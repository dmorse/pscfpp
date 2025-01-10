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
#include <pscf/cuda/GpuResources.h>
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
      if (!isAllocated_) {
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
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
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b); 
   }

   // Compute and return inner product of two vectors.
   template <int D>
   double AmCompressor<D>::dotProduct(DeviceArray<cudaReal> const & a, 
                                      DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double AmCompressor<D>::maxAbs(DeviceArray<cudaReal> const & a)
   {
      return Reduce::maxAbs(a);
   }

   // Update basis
   template <int D>
   void 
   AmCompressor<D>::updateBasis(RingBuffer< DeviceArray<cudaReal> > & basis,
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
   AmCompressor<D>::addHistories(DeviceArray<cudaReal>& trial,
                                 RingBuffer<DeviceArray<cudaReal> > const & basis,
                                 DArray<double> coeffs,
                                 int nHist)
   {
      for (int i = 0; i < nHist; i++) {
         VecOp::addEqVc(trial, basis[i], -1.0 * coeffs[i]);
      }
   }

   template <int D>
   void AmCompressor<D>::addPredictedError(DeviceArray<cudaReal>& fieldTrial,
                                           DeviceArray<cudaReal> const & resTrial,
                                           double lambda)
   {
      VecOp::addEqVc(fieldTrial, resTrial, lambda);
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
      
      // Initialize resid to c field of species 0 minus 1
      VecOp::subVS(resid, system().c().rgrid(0), 1.0);

      // Add other c fields to get SCF residual vector elements
      for (int i = 1; i < nMonomer; i++) {
         VecOp::addEqV(resid, system().c().rgrid(i));
      }
   }

   // Update the current system field coordinates
   template <int D>
   void AmCompressor<D>::update(DeviceArray<cudaReal>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      
      // New field is the w0_ + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(wFieldTmp_[i], w0_[i], newGuess);
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