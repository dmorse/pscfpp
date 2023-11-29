#ifndef PSPC_LR_AM_COMPRESSOR_TPP
#define PSPC_LR_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <pscf/chem/Monomer.h>
#include <pscf/mesh/MeshIterator.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <pspc/System.h>
#include <util/global.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   // Constructor
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
    : Compressor<D>(system),
      isAllocated_(false)
   {  setClassName("LrAmCompressor"); }

   // Destructor
   template <int D>
   LrAmCompressor<D>::~LrAmCompressor()
   {  setClassName("LrAmCompressor"); }

   // Read parameters from file
   template <int D>
   void LrAmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readErrorType(in);
   }
      
   // Initialize just before entry to iterative loop.
   template <int D>
   void LrAmCompressor<D>::setup(bool isContinuation)
   {  
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DArray<double> >::setup(isContinuation);
      
      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
         }
      }
      
      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_){
         newBasis_.allocate(meshSize);
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         error_.allocate(meshSize);
         resid_.allocate(dimensions);
         residK_.allocate(dimensions);
         intraCorrelation_.allocate(kMeshDimensions_);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(meshSize);
            wFieldTmp_[i].allocate(meshSize);
         }
            
         isAllocated_ = true;
      }
      
      // Store value of initial guess chemical potential fields
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j< meshSize; ++j){
            w0_[i][j] = system().w().rgrid(i)[j];
         }
      }
      
      computeIntraCorrelation();
   }
  
   // Iterative solver (AM algorithm) 
   template <int D>
   int LrAmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DArray<double> >::solve();
      mdeCounter_ = AmIteratorTmpl<Compressor<D>,DArray<double>>::totalItr();
      return solve;
   }

   // Assign one array to another
   template <int D>
   void LrAmCompressor<D>::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

   // Compute and return inner product of two vectors.
   template <int D>
   double LrAmCompressor<D>::dotProduct(DArray<double> const & a, 
                                    DArray<double> const & b)
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = 0.0;
      for (int i = 0; i < n; i++) {
         // If either value is NaN, throw NanException
         if (std::isnan(a[i]) || std::isnan(b[i])) { 
            throw NanException("LrAmCompressor::dotProduct",
                               __FILE__,__LINE__,0);
         }
         product += a[i] * b[i];
      }
      return product;
   }

   // Compute and return maximum element of a vector.
   template <int D>
   double LrAmCompressor<D>::maxAbs(DArray<double> const & a)
   {
      const int n = a.capacity();
      double max = 0.0;
      double value;
      for (int i = 0; i < n; i++) {
         value = a[i];
         if (std::isnan(value)) { // if value is NaN, throw NanException
            throw NanException("LrAmCompressor::dotProduct",__FILE__,__LINE__,0);
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
   LrAmCompressor<D>::updateBasis(RingBuffer< DArray<double> > & basis,
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
   LrAmCompressor<D>::addHistories(DArray<double>& trial,
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
   void 
   LrAmCompressor<D>::addPredictedError(DArray<double>& fieldTrial,
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
   bool LrAmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return the number of elements in a field vector
   template <int D>
   int LrAmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   // Get the current field from the system
   template <int D>
   void LrAmCompressor<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into  linear arrays
      const int meshSize = system().domain().mesh().size();
      const DArray< RField<D> > * currSys = &system().w().rgrid();
      
      /*
      * The field that we are adjusting is the Langrange multiplier field 
      * defined on a grid.  The current value is the difference between 
      * w and w0_ for the first monomer (or any other monomer).
      */
      for (int i = 0; i < meshSize; i++){
         curr[i] = (*currSys)[0][i] - w0_[0][i];
      }

   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void LrAmCompressor<D>::evaluate()
   {  system().compute(); }

   // Compute the residual for the current system state
   template <int D>
   void LrAmCompressor<D>::getResidual(DArray<double>& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const double vMonomer = system().mixture().vMonomer();  
      
      // Initialize residuals
      for (int i = 0 ; i < n; ++i) {
         resid_[i] = -1.0;
      }

      // Compute SCF residual vector elements
      for (int j = 0; j < nMonomer; ++j) {
        for (int k = 0; k < meshSize; ++k) {
           resid_[k] += system().c().rgrid(j)[k];
        }
      }
      
      // Convert residual to Fourier Space
      system().fft().forwardTransform(resid_, residK_);
      // Residual combine with Linear response factor
      MeshIterator<D> iter;
      iter.setDimensions(residK_.dftDimensions());
      for (iter.begin(); !iter.atEnd(); ++iter) {
         residK_[iter.rank()][0] *= 1.0 / (vMonomer * intraCorrelation_[iter.rank()]);
         residK_[iter.rank()][1] *= 1.0 / (vMonomer * intraCorrelation_[iter.rank()]);
      }
      
      // Convert back to real Space
      system().fft().inverseTransform(residK_, resid_);
      
      // Copy to resid
      for (int i = 0 ; i < n; ++i) {
         resid[i] = resid_[i];
      }
   }

   // Update the current system field coordinates
   template <int D>
   void LrAmCompressor<D>::update(DArray<double>& newGuess)
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
   void LrAmCompressor<D>::outputToLog()
   {}
   
   template<int D>
   void LrAmCompressor<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Compressor times contributions:\n";
      AmIteratorTmpl<Compressor<D>, DArray<double> >::outputTimers(out);
   }
   
   
   template<int D>
   void LrAmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DArray<double> >::clearTimers();
      mdeCounter_ = 0;
   }
   
   template<int D>
   double LrAmCompressor<D>::computeDebye(double x)
   {
      if (x == 0){
         return 1.0;
      } else {
         return 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
   }
   
   template<int D>
   double LrAmCompressor<D>::computeIntraCorrelation(double qSquare)
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
         for (int j = 0; j < nBlock; j++) {
            monomerId = polymerPtr-> block(j).monomerId();
            kuhn = system().mixture().monomer(monomerId).kuhn();
            // Get the length (number of monomers) in this block.
            length = polymerPtr-> block(j).length();
            rg2 = length * kuhn* kuhn /6.0;
            g = computeDebye(qSquare*rg2);
            omega += length * g/ vMonomer;
         }
      }
      return omega;
   }
   
   template<int D>
   void LrAmCompressor<D>::computeIntraCorrelation()
   {
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, system().mesh().dimensions(), 
                                  system().unitCell());
         Gsq = system().unitCell().ksq(Gmin);
         intraCorrelation_[iter.rank()] = computeIntraCorrelation(Gsq);
      }
   }
   
   template<int D>
   double LrAmCompressor<D>::computeError(int verbose)
   {
      errorType_ = AmIteratorTmpl<Compressor<D>, DArray<double> >::errorType();
      double error = 0.0;
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      // Initialize residuals
      for (int i = 0 ; i < n; ++i) {
         error_[i] = -1.0;
      }

       // Compute SCF residual vector elements
      for (int j = 0; j < nMonomer; ++j) {
        for (int k = 0; k < meshSize; ++k) {
         error_[k] += system().c().rgrid(j)[k];
        }
      }
      
      // Find max residual vector element
      double maxRes  = maxAbs(error_);
     
      // Find norm of residual vector
      double normRes = AmIteratorTmpl<Compressor<D>, DArray<double> >::norm(error_);
      
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
         if (errorType_ == "maxResid") {
            error = maxRes;
         } else if (errorType_ == "normResid") {
            error = normRes;
         } else if (errorType_ == "rmsResid") {
            error = rmsRes;
         } else {
            UTIL_THROW("Invalid iterator error type in parameter file.");
         }

      } else {
         // Set error value
         if (errorType_ == "maxResid") {
            error = maxRes;
         } else if (errorType_ == "normResid") {
            error = normRes;
         } else if (errorType_ == "rmsResid") {
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
