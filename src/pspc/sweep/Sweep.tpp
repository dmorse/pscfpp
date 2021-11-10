/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include <pspc/System.h>
#include <pspc/iterator/AmIterator.h>
#include <pscf/sweep/SweepTmpl.tpp>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   Sweep<D>::Sweep() 
    : SweepTmpl< BasisFieldState<D> >(2),
      systemPtr_(0)
   {}

   /*
   * Constructor, creates association with parent system.
   */
   template <int D>
   Sweep<D>::Sweep(System<D> & system) 
    : SweepTmpl< BasisFieldState<D> >(2),
      systemPtr_(0)
   {  setSystem(system); }

   /*
   * Destructor.
   */
   template <int D>
   Sweep<D>::~Sweep() 
   {}

   /*
   * Set association with a parent system.
   */
   template <int D>
   void Sweep<D>::setSystem(System<D> & system) 
   {  systemPtr_ = & system; }

   /*
   * Check allocation of one state object, allocate if necessary.
   */
   template <int D>
   void Sweep<D>::checkAllocation(BasisFieldState<D>& state) 
   {  state.allocate(); }

   /*
   * Setup operations at the beginning of a sweep.
   */
   template <int D>
   void Sweep<D>::setup() 
   {
      initialize();

      // Check or create associations of states with the parent System
      // Note: FieldState::setSystem does nothing if already set.
      for (int i=0; i < historyCapacity(); ++i) {
         state(i).setSystem(system());
      }
      trial_.setSystem(system());
   };

   /*
   * Set non-adjustable system parameters to new values.
   *
   * \param s path length coordinate, in range [0,1]
   */
   template <int D>
   void Sweep<D>::setParameters(double s) 
   {
      // Empty default implementation allows Sweep<D> to be compiled.
      UTIL_THROW("Calling unimplemented function Sweep::setParameters");
   };

   /*
   * Create guess for adjustable variables by continuation.
   */
   template <int D>
   void Sweep<D>::setGuess(double sNew) 
   {

      // If historySize() == 1, do nothing, use previous system state
      // as guess for the new state.

      if (historySize() > 1) {

         UTIL_CHECK(historySize() <= historyCapacity());

         // Compute coefficients of previous states for continuation
         setCoefficients(sNew);

         bool isFlexible = system().iterator().isFlexible();

         // Set trial w fields 
         int nMonomer = system().mixture().nMonomer();
         int nStar = system().basis().nStar();
         DArray<double>* newFieldPtr;
         DArray<double>* oldFieldPtr; 
         int i, j, k;
         for (i=0; i < nMonomer; ++i) {
            newFieldPtr = &trial_.field(i);

            // Previous state k = 0 (most recent)
            oldFieldPtr = &state(0).field(i);
            for (j=0; j < nStar; ++j) {
               (*newFieldPtr)[j] = (*oldFieldPtr)[j]*c(0);
            }

            // Previous states k >= 1 (older)
            for (k = 1; k < historySize(); ++k) {
               oldFieldPtr = &state(k).field(i);
               for (j=0; j < nStar; ++j) {
                  (*newFieldPtr)[j] += (*oldFieldPtr)[j]*c(k);
               }
            }
         }

         // Set trial unit cell
         if (isFlexible) {
            int nParameter = system().unitCell().nParameter();
            unitCellParameters_.clear();
            for (int i = 0; i < nParameter; ++i) {
               unitCellParameters_.append(state(0).unitCell().parameter(i)*c(0));
            }
            for (k = 1; k < historySize(); ++k) {
               for (int i = 0; i < nParameter; ++i) {
                  unitCellParameters_[i] += state(0).unitCell().parameter(i)*c(k);
               }
            }
            trial_.unitCell().setParameters(unitCellParameters_);
         }

         trial_.setSystemState(isFlexible);
      }

    
   };

   /*
   * Call current iterator to solve SCFT problem.
   *
   * Return 0 for sucessful solution, 1 on failure to converge.
   */
   template <int D>
   int Sweep<D>::solve(bool isContinuation)
   {  return system().iterate(); };

   /*
   * Reset system to previous solution after iterature failure.
   *
   * The implementation of this function should reset the system state
   * to correspond to that stored in state(0).
   */
   template <int D>
   void Sweep<D>::reset()
   {  
      bool isFlexible = system().iterator().isFlexible();
      state(0).setSystemState(isFlexible); 
   }

   /**
   * Update state(0) and output data after successful convergence
   *
   * The implementation of this function should copy the current 
   * system state into state(0) and output any desired information
   * about the current converged solution.
   */
   template <int D>
   void Sweep<D>::getSolution() 
   { 
      state(0).getSystemState(); 
      // Output operations (see example in Fd1d)
   };

} // namespace Pspc
} // namespace Pscf
