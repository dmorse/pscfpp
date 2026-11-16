#ifndef RP_POLYMER_TPP
#define RP_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>
#include <pscf/field/RField.h>
#include <pscf/backend/Reduce.h>
#include <pscf/backend/VecOp.h>
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>

#include <cmath>
#include <vector>

namespace Pscf {
namespace Rp {

   /*
   * Constructor.
   */
   template <int D, class T>
   Polymer<D,T>::Polymer()
    : stress_(),
      nParam_(0)
   {  ParamComposite::setClassName("Polymer"); }

   /*
   * Destructor.
   */
   template <int D, class T>
   Polymer<D,T>::~Polymer()
   {}

   /*
   * Set the number of unit cell parameters.
   */
   template <int D, class T>
   void Polymer<D,T>::setNParams(int nParam)
   {  nParam_ = nParam; }

   /*
   * Clear all data that depends on unit cell parameters.
   */
   template <int D, class T>
   void Polymer<D,T>::clearUnitCellData()
   {
      for (int j = 0; j < nBlock(); ++j) {
         block(j).clearUnitCellData();
      }
      stress_.clear();
   }

   /*
   * Compute solution to MDE and block concentrations.
   */
   template <int D, class T>
   void Polymer<D,T>::compute(DArray< RField<D,T> > const & wFields,
                              double phiTot)
   {
      // Setup solvers for all blocks. Each block solves the MDE in the
      // field w of its monomer type minus the shift
      // min(w) - ln<exp(-n(w - min(w)))>/n (see Polymer.h), which
      // PolymerTmpl::solve removes from ln q analytically.
      const double n = PolymerModel::isThread() ? length()
                                                : (double) nBead();
      const int nb = nBlock();
      std::vector<double> shifts(nb);
      for (int j = 0; j < nb; ++j) {
         RField<D,T> const & w = wFields[block(j).monomerId()];
         if (wWork_.capacity() != w.capacity()) {
            if (wWork_.isAllocated()) wWork_.deallocate();
            wWork_.allocate(w.meshDimensions());
         }
         const double wMin = Reduce::min(w);
         VecOp::subVS(wWork_, w, wMin);
         VecOp::expVc(wWork_, wWork_, -n);
         shifts[j] = wMin
                  - std::log(Reduce::sum(wWork_)/double(w.capacity()))/n;
         block(j).setupSolver(w, shifts[j]);
      }

      // Call base class PolymerTmpl solve() function
      // This solves the MDE for all propagators in a precalculated order
      PolymerTmplT::solve(phiTot);

      // If qShifted() is too small, raise the shifts (see raiseShifts)
      auto setupBlock = [&](int j, double s) {
         block(j).setupSolver(wFields[block(j).monomerId()], s);
      };
      auto solveBlock = [&](int j) -> double {
         VecOp::eqS(wWork_, 1.0);
         PropagatorT& p = block(j).propagator(0);
         p.solve(wWork_);
         // The last slice (tail() is not defined for a bead chain end)
         return Reduce::sum(p.q(p.ns()-1))/double(wWork_.capacity());
      };
      PolymerTmplT::raiseShifts(shifts, setupBlock, solveBlock, phiTot);

      // Compute block concentration fields, normalizing products of
      // the shifted propagators by their own partition function
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( qShifted() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationThread(prefactor);
         }
      } else {
         prefactor = phi() / ( qShifted() * (double)nBead() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationBead(prefactor);
         }
      }

   }

   /*
   * Compute stress contribution from a polymer species.
   */
   template <int D, class T>
   void Polymer<D,T>::computeStress()
   {
      UTIL_CHECK(nParam_ > 0);

      // Compute stress contributions for all blocks
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( qShifted() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeStressThread(prefactor);
         }
      } else {
         prefactor = phi() / ( qShifted() * (double)nBead() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeStressBead(prefactor);
         }
      }

      // Initialize all stress_ elements to zero
      stress_.clear();
      for (int i = 0; i < nParam_; ++i) {
         stress_.append(0.0);
      }

      // Sum over all block stress contributions
      for (int i = 0; i < nBlock(); ++i) {
         for (int j = 0; j < nParam_; ++j){
            stress_[j] += block(i).stress(j);
         }
      }

   }

}
}
#endif
