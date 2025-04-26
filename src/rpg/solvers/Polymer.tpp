#ifndef RPG_POLYMER_TPP
#define RPG_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rpg/solvers/Block.h>
#include <prdc/cuda/RField.h>
#include <pscf/chem/PolymerModel.h>

namespace Pscf {
namespace Rpg { 

   /*
   * Constructor.
   */
   template <int D>
   Polymer<D>::Polymer()
    : nParams_(0)
   {  ParamComposite::setClassName("Polymer"); }

   /*
   * Destructor.
   */
   template <int D>
   Polymer<D>::~Polymer()
   {}

   /*
   * Store the number of lattice parameters in the unit cell.
   */ 
   template <int D>
   void Polymer<D>::setNParams(int nParams)
   {  nParams_ = nParams; }

   /*
   * Compute solution to MDE and concentrations.
   */ 
   template <int D>
   void Polymer<D>::compute(DArray< RField<D> > const & wFields, 
                            double phiTot)
   {
      // Setup solvers for all blocks
      int monomerId;
      for (int j = 0; j < nBlock(); ++j) {
         monomerId = block(j).monomerId();
         block(j).setupSolver(wFields[monomerId]);
      }

      // Call base class PolymerTmpl solve() function
      // This solves the MDE for all propagators in a precalculated order
      solve(phiTot);

      // Compute block concentration fields
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( q() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationThread(prefactor);
         }
      } else {
         prefactor = phi() / ( q() * (double)nBead() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationBead(prefactor);
         }
      }

   }

   /*
   * Compute stress contribution from a polymer species.
   */
   template <int D>
   void Polymer<D>::computeStress()
   {
      UTIL_CHECK(nParams_ > 0);
      
      // Compute stress contributions for all blocks
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( q() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeStressThread(prefactor);
         }
      } else {
         prefactor = phi() / ( q() * (double)nBead() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeStressBead(prefactor);
         }
      }

      // Initialize all stress_ elements to zero
      for (int i = 0; i < nParams_; ++i) {
        stress_[i] = 0.0;
      }

      // Sum over all block stress contributions
      for (int i = 0; i < nBlock(); ++i) {
         for (int j=0; j < nParams_; ++j){
            stress_[j] += block(i).stress(j);
         }
      }

   }

}
}
#endif
