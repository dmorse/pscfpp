#ifndef RPC_AM_COMPRESSOR_TPP
#define RPC_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
    : Compressor<D>(system),
      isAllocated_(false)
   {  ParamComposite::setClassName("AmCompressor"); }

   /*
   * Destructor.
   */
   template <int D>
   AmCompressor<D>::~AmCompressor()
   {}

   /*
   * Read parameters from file.
   */
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      using AmTmpl = AmIteratorTmpl<Compressor<D>, DArray<double> >;
      bool useLambdaRamp = false;  // Default value

      AmTmpl::readParameters(in);
      AmTmpl::readErrorType(in);
      AmTmpl::readMixingParameters(in, useLambdaRamp);
   }

   /*
   * Initialize just before entry to iterative loop.
   */
   template <int D>
   void AmCompressor<D>::setup(bool isContinuation)
   {
      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>, DArray<double> >::setup(isContinuation);

      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();

      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_){
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

   /*
   * Main function - identify partial saddle-point state.
   */
   template <int D>
   int AmCompressor<D>::compress()
   {
      int solve = AmIteratorTmpl<Compressor<D>, DArray<double> >::solve();
      return solve;
   }

   /*
   * Output timer information, if requested.
   */
   template<int D>
   void AmCompressor<D>::outputTimers(std::ostream& out) const
   {
      // Output timing results, if requested.
      out << "\n";
      out << "AmCompressor time contributions:\n";
      AmIteratorTmpl<Compressor<D>, DArray<double> >::outputTimers(out);
   }

   /*
   * Clear timers and MDE counter.
   */
   template<int D>
   void AmCompressor<D>::clearTimers()
   {
      AmIteratorTmpl<Compressor<D>, DArray<double> >::clearTimers();
      mdeCounter_ = 0;
   }

   // Private virtual functions that interact with parent system

   /*
   * Does the system have an initial field guess?
   */
   template <int D>
   bool AmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D>
   int AmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Get the current field from the system.
   */
   template <int D>
   void AmCompressor<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into  linear arrays
      const int meshSize = system().domain().mesh().size();
      const DArray< RField<D> > * currSys = &system().w().rgrid();

      /*
      * The field that we are adjusting is the Langrange multiplier 
      * field.  The current value is the difference between w and w0_ 
      * for the first monomer type, but any monomer type would give
      * the same answer.
      */
      for (int i = 0; i < meshSize; i++){
         curr[i] = (*currSys)[0][i] - w0_[0][i];
      }

   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void AmCompressor<D>::evaluate()
   {
      system().compute();
      ++mdeCounter_;
   }

   /*
   * Compute the residual vector for the current system state.
   */
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

   /*
   * Update the current system field coordinates.
   */
   template <int D>
   void AmCompressor<D>::update(DArray<double>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();

      // New field is the w0_ + newGuess for the pressure field
      for (int i = 0; i < nMonomer; i++){
         for (int k = 0; k < meshSize; k++){
            wFieldTmp_[i][k] = w0_[i][k] + newGuess[k];
         }
      }
      system().w().setRGrid(wFieldTmp_);
   }

   /*
   * Do-nothing output function.
   */
   template<int D>
   void AmCompressor<D>::outputToLog()
   {}

}
}
#endif
