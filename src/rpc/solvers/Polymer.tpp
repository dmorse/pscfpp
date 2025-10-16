#ifndef RPC_POLYMER_TPP
#define RPC_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Block.h"
#include "Propagator.h"
#include <pscf/solvers/PolymerTmpl.tpp>
#include <prdc/cpu/RField.h>

namespace Pscf {
namespace Rpc { 

   /*
   * Constructor.
   */
   template <int D>
   Polymer<D>::Polymer()
    : stress_(),
      nParam_(0)
   {  ParamComposite::setClassName("Polymer");}

   /*
   * Destructor.
   */
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
   * Clear all data that depends on unit cell parameters.
   */ 
   template <int D>
   void Polymer<D>::clearUnitCellData()
   {
      for (int j = 0; j < nBlock(); ++j) {
         block(j).clearUnitCellData();
      }
      stress_.clear();
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
      stress_.clear();
      for (int i = 0; i < nParam_; ++i) {
        stress_.append(0.0);
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
