/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Block.h"
#include "Propagator.h"
#include <r1/field/Domain.h>
#include <pscf/solvers/PolymerTmpl.tpp>

#include <cmath>
#include <vector>

// Explicit instantiation of base class
namespace Pscf {
   template class PolymerTmpl<R1::Block, R1::Propagator, double>;
}

namespace Pscf {
namespace R1 { 

   /*
   * Constructor
   */
   Polymer::Polymer()
   {  setClassName("Polymer"); }

   /*
   * Destructor
   */
   Polymer::~Polymer()
   {}

   /*
   * Compute solution to MDE and concentrations.
   */ 
   void Polymer::compute(DArray< DArray<double> > const & wFields)
   {

      // Setup solvers for all blocks. Each block solves the MDE in the
      // field w of its monomer type minus the shift
      // min(w) - ln<exp(-N(w - min(w)))>/N, N = length() (see Polymer.h),
      // which PolymerTmpl::solve removes from ln q analytically.
      const double n = length();
      const int nb = nBlock();
      std::vector<double> shifts(nb);
      for (int j = 0; j < nb; ++j) {
         DArray<double> const & w = wFields[block(j).monomerId()];
         const int nx = w.capacity();
         if (wWork_.capacity() != nx) {
            if (wWork_.isAllocated()) wWork_.deallocate();
            wWork_.allocate(nx);
         }
         double wMin = w[0];
         for (int i = 1; i < nx; ++i) {
            if (w[i] < wMin) wMin = w[i];
         }
         for (int i = 0; i < nx; ++i) {
            wWork_[i] = std::exp(-n*(w[i] - wMin));
         }
         shifts[j] = wMin
                    - std::log(block(j).domain().spatialAverage(wWork_))/n;
         block(j).setupSolver(w, shifts[j]);
      }

      // Solve MDE for all propagators
      solve();

      // If qShifted() is too small, raise the shifts (see raiseShifts)
      auto setupBlock = [&](int j, double s) {
         block(j).setupSolver(wFields[block(j).monomerId()], s);
      };
      auto solveBlock = [&](int j) -> double {
         for (int i = 0; i < wWork_.capacity(); ++i) {
            wWork_[i] = 1.0;
         }
         Propagator& p = block(j).propagator(0);
         p.solve(wWork_);
         return block(j).domain().spatialAverage(p.tail());
      };
      raiseShifts(shifts, setupBlock, solveBlock);

      // Compute block concentration fields (thread model), normalizing
      // products of the shifted propagators by their own partition
      // function
      double prefactor;
      prefactor = phi() / ( qShifted() * length() );
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }

   }

}
}
