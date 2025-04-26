#ifndef RPC_AM_COMPRESSOR_TPP
#define RPC_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpc/System.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc{

   using namespace Util;

   // Constructor
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
    : Compressor<D>(system),
      isAllocated_(false)
   {  setClassName("AmCompressor"); }

   // Destructor
   template <int D>
   AmCompressor<D>::~AmCompressor()
   {}

   // Read parameters from file
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readErrorType(in);
   }

   // Initialize just before entry to iterative loop.
   template <int D>
   void AmCompressor<D>::setup(bool isContinuation)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DArray<double> >::setup(isContinuation);

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
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j< meshSize; ++j){
            w0_[i][j] = system().w().rgrid(i)[j];
         }
      }
   }

   template <int D>
   int AmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DArray<double> >::solve();
      //mdeCounter_ = AmIteratorTmpl<Compressor<D>,DArray<double>>::totalItr();
      return solve;
   }

   // Assign one array to another
   template <int D>
   void AmCompressor<D>::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

   // Compute and return inner product of two vectors.
   template <int D>
   double AmCompressor<D>::dotProduct(DArray<double> const & a,
                                    DArray<double> const & b)
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

   // Compute and return maximum element of a vector.
   template <int D>
   double AmCompressor<D>::maxAbs(DArray<double> const & a)
   {
      const int n = a.capacity();
      double max = 0.0;
      double value;
      for (int i = 0; i < n; i++) {
         value = a[i];
         if (std::isnan(value)) { // if value is NaN, throw NanException
            throw NanException("AmCompressor::dotProduct",__FILE__,__LINE__,0);
         }
         if (fabs(value) > max) {
            max = fabs(value);
         }
      }
      return max;
   }

   // Update basis
   template <int D>
   void
   AmCompressor<D>::updateBasis(RingBuffer< DArray<double> > & basis,
                              RingBuffer< DArray<double> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();

      // New basis vector is difference between two most recent states
      for (int i = 0; i < n; i++) {
         newBasis_[i] = hists[0][i] - hists[1][i];
      }

      basis.append(newBasis_);
   }

   template <int D>
   void
   AmCompressor<D>::addHistories(DArray<double>& trial,
                               RingBuffer<DArray<double> > const & basis,
                               DArray<double> coeffs,
                               int nHist)
   {
      int n = trial.capacity();
      for (int i = 0; i < nHist; i++) {
         for (int j = 0; j < n; j++) {
            // Not clear on the origin of the -1 factor
            trial[j] += coeffs[i] * -1 * basis[i][j];
         }
      }
   }

   template <int D>
   void AmCompressor<D>::addPredictedError(DArray<double>& fieldTrial,
                                         DArray<double> const & resTrial,
                                         double lambda)
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
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
   void AmCompressor<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into  linear arrays
      const int meshSize = system().domain().mesh().size();
      const DArray< RField<D> > * currSys = &system().w().rgrid();

      /*
      * The field that we are adjusting is the Langrange multiplier field
      * with number of grid pts components.The current value is the difference
      * between w and w0_ for the first monomer (any monomer should give the same answer)
      */
      for (int i = 0; i < meshSize; i++){
         curr[i] = (*currSys)[0][i] - w0_[0][i];
      }

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
   void AmCompressor<D>::getResidual(DArray<double>& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();

      // Initialize residuals
      for (int i = 0 ; i < n; ++i) {
         resid[i] = -1.0;
      }

       // Compute SCF residual vector elements
      for (int j = 0; j < nMonomer; ++j) {
        for (int k = 0; k < meshSize; ++k) {
           resid[k] += system().c().rgrid(j)[k];
        }
      }

   }

   // Update the current system field coordinates
   template <int D>
   void AmCompressor<D>::update(DArray<double>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();

      //New field is the w0_ + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            wFieldTmp_[i][k] = w0_[i][k] + newGuess[k];
         }
      }
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
      out << "AmCompressor time contributions:\n";
      AmIteratorTmpl<Compressor<D>, DArray<double> >::outputTimers(out);
   }

   // Clear timers and MDE counter
   template<int D>
   void AmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DArray<double> >::clearTimers();
      mdeCounter_ = 0;
   }

}
}
#endif
