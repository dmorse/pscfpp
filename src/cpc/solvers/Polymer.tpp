#ifndef CPC_POLYMER_TPP
#define CPC_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Block.h"
#include "Propagator.h"
#include <prdc/field/cpu/CField.h>
#include <pscf/cpu/complex.h>

namespace Pscf {
namespace Cpc { 

   /*
   * Constructor.
   */
   template <int D>
   Polymer<D>::Polymer()
   {  ParamComposite::setClassName("Polymer");}

   /*
   * Destructor.
   */
   template <int D>
   Polymer<D>::~Polymer()
   {}

   /*
   * Clear all data that depends on unit cell parameters.
   */ 
   template <int D>
   void Polymer<D>::clearUnitCellData()
   {
      for (int j = 0; j < nBlock(); ++j) {
         block(j).clearUnitCellData();
      }
   }

   /*
   * Compute solution to MDE and block concentrations.
   */ 
   template <int D>
   void Polymer<D>::compute(DArray< CField<D> > const & wFields)
   {
      // Setup solvers for all blocks
      int monomerId;
      for (int j = 0; j < nBlock(); ++j) {
         monomerId = block(j).monomerId();
         block(j).setupSolver(wFields[monomerId]);
      }

      // Call base class PolymerTmpl solve() function
      // Solve MDE for all propagators
      double phiTot = 1.0;
      solve(phiTot);

      // Compute block concentration fields
      std::complex<double> ratio;
      fftw_complex prefactor;
      if (PolymerModel::isThread()) {
         ratio = phi() / length();
         ratio /= q();
         assign(prefactor, ratio);
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationThread(prefactor);
         }
      } else {
         UTIL_CHECK(PolymerModel::isBead());
         double len = (double) nBead();
         ratio = phi() / len;
         ratio /= q();
         assign(prefactor, ratio);
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationBead(prefactor);
         }
      }

   }

}
}
#endif
