#ifndef RPC_LR_AM_COMPRESSOR_TPP
#define RPC_LR_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LrAmCompressor.h"
#include <rpc/system/System.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc{

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   LrAmCompressor<D>::LrAmCompressor(System<D>& system)
    : Compressor<D>(system),
      intra_(system),
      isIntraCalculated_(false),
      isAllocated_(false)
   {  setClassName("LrAmCompressor"); }

   /*
   * Destructor.
   */
   template <int D>
   LrAmCompressor<D>::~LrAmCompressor()
   {}

   /*
   * Read parameters from file.
   */
   template <int D>
   void LrAmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readErrorType(in);
   
   }

   /*
   * Initialize just before entry to iterative loop.
   */
   template <int D>
   void LrAmCompressor<D>::setup(bool isContinuation)
   {
      
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DArray<double> >::setup(isContinuation);
      
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      
      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
         }
      }
      
      // Allocate memory required by compressor, if not done previously
      if (!isAllocated_){
         newBasis_.allocate(meshSize);
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
      
      // Store initial values of monomer chemical potential fields in w0_
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j< meshSize; ++j){
            w0_[i][j] = system().w().rgrid(i)[j];
         }
      }
      
   }

   /*
   * Apply the Anderson-Mixing algorithm (main function).
   */
   template <int D>
   int LrAmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DArray<double> >::solve();
      return solve;
   }

   /*
   * Output timing results to an output stream (file)
   */
   template<int D>
   void LrAmCompressor<D>::outputTimers(std::ostream& out) const
   {
      // Output timing results, if requested.
      out << "\n";
      out << "LrAmCompressor time contributions:\n";
      AmIteratorTmpl<Compressor<D>, DArray<double> >::outputTimers(out);
   }

   /*
   * Clear all timers and the MDE counter
   */
   template<int D>
   void LrAmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DArray<double> >::clearTimers();
      mdeCounter_ = 0;
   }


   // Private virtual function that interact with parent system

   /*
   * Correction step (second step of Anderson mixing)
   *
   * This LrAM algorithm uses a quasi-Newton correction step with an
   * approximate Jacobian given by the Jacobian in a homogeneous state.
   */
   template <int D>
   void
   LrAmCompressor<D>::addPredictedError(DArray<double>& fieldTrial,
                                        DArray<double> const & resTrial)
   {
      // Local constants
      const int n = fieldTrial.capacity();
      const double vMonomer = system().mixture().vMonomer();

      // Copy DArray<double> resTrial into RField<D> resid_ 
      // Allows use of FFT functions that take an RField container
      for (int i = 0 ; i < n; ++i) {
         resid_[i] = resTrial[i];
      }

      // Fourier transform r-grid residual resid_ to k-grid form residK_
      system().domain().fft().forwardTransform(resid_, residK_);

      // Compute update on a k-grid, using quasi-Newton algorithm
      double factor;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         factor = 1.0 / (vMonomer * intraCorrelationK_[iter.rank()]);
         residK_[iter.rank()][0] *= factor;
         residK_[iter.rank()][1] *= factor;
      }
      // Field residK_ now stores k-grid update rather than residual

      // Convert field update back to r-grid (real space) form 
      // On return, resid_ contains r-grid field update, residK_ is destroyed
      system().domain().fft().inverseTransformUnsafe(residK_, resid_);

      // Add update to obtain new fieldTrial
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += resid_[i];
      }
   }

   // Private virtual functions that interact with parent System

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D>
   int LrAmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Does the associated system have initialized w fields?
   */
   template <int D>
   bool LrAmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field variable from the system.
   *
   * The field variable is the change in the Lagrange multiplier field.
   * relative to that used in initial array of monomer fields, w0_.
   */
   template <int D>
   void LrAmCompressor<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into  linear arrays
      const int meshSize = system().domain().mesh().size();
      const DArray< RField<D> > * currSys = &system().w().rgrid();

      /*
      * The field that we are adjusting is the Langrange multiplier field
      * with number of grid pts components. The current value is the
      * difference between w and w0_ for the first monomer (any monomer
      * should give the same answer)
      */
      for (int i = 0; i < meshSize; i++){
         curr[i] = (*currSys)[0][i] - w0_[0][i];
      }

   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void LrAmCompressor<D>::evaluate()
   {
      system().compute();
      ++mdeCounter_;
   }

   /*
   * Compute the residual vector for the current system state
   */
   template <int D>
   void LrAmCompressor<D>::getResidual(DArray<double>& resid)
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

   /*
   * Update the w field values stored in the system 
   */
   template <int D>
   void LrAmCompressor<D>::update(DArray<double>& newGuess)
   {
      // Local constant values
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();

      // Locally update monomer w fields in array wFieldTmp_
      // Set wFieldTmp = w0_ + newGuess 
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            wFieldTmp_[i][k] = w0_[i][k] + newGuess[k];
         }
      }

      // Update w fields in the system WContainer
      system().w().setRGrid(wFieldTmp_);
   }

   /*
   * Output results to log file (do-nothing implementation).
   */
   template<int D>
   void LrAmCompressor<D>::outputToLog()
   {}

}
}
#endif
