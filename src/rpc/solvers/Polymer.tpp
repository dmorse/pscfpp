#ifndef RPC_POLYMER_TPP
#define RPC_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pscf {
namespace Rpc { 

   template <int D>
   Polymer<D>::Polymer()
    : stress_(),
      nParam_(0)
   {  ParamComposite::setClassName("Polymer");}

   template <int D>
   Polymer<D>::~Polymer()
   {}

   /*
   * Set the number of unit cell parameters.
   */ 
   template <int D>
   void Polymer<D>::setNParams(int nParam)
   {   nParam_ = nParam; }

   /*
   * Set unit cell dimensions in all solvers.
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
      // Solve MDE for all propagators
      solve(phiTot);

      // Compute block concentration fields
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi() / ( q() * length() );
         for (int i = 0; i < nBlock(); ++i) {
            block(i).computeConcentrationThread(prefactor);
         }
      } else 
      if (PolymerModel::isBead()) {
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

      // Initialize all stress_ elements zero
      for (int i = 0; i < 6; ++i) {
        stress_[i] = 0.0;
      }

      // Sum over all block stress contributions
      for (int i = 0; i < nBlock(); ++i) {
         for (int j = 0; j < nParam_ ; ++j){
            stress_[j] += block(i).stress(j);
         }
      }

   }

}
}
#endif
