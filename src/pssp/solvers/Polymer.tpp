#ifndef PSSP_POLYMER_TPP
#define PSSP_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pscf {
namespace Pssp { 

   template <int D>
   Polymer<D>::Polymer()
   {
      setClassName("Polymer"); 
   }

   template <int D>
   Polymer<D>::~Polymer()
   {}

   template <int D>
   void Polymer<D>::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Species::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi; 
   }

   template <int D>
   void Polymer<D>::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Species::Open);  
      mu_ = mu; 
   }

   /*
   * Set unit cell dimensions in all solvers.
   */ 
   template <int D>
   void Polymer<D>::setupUnitCell(UnitCell<D> const & unitCell)
   {
      // Set association to unitCell
      unitCellPtr_ = &unitCell;

      for (int j = 0; j < nBlock(); ++j) {
         block(j).setupUnitCell(unitCell);
      }
   }

   /*
   * Compute solution to MDE and concentrations.
   */ 
   template <int D>
   void Polymer<D>::compute(DArray<WField> const & wFields)
   {
      // Setup solvers for all blocks
      int monomerId;
      for (int j = 0; j < nBlock(); ++j) {
         monomerId = block(j).monomerId();
         block(j).setupSolver(wFields[monomerId]);
      }

      // Call generic solver() method base class template.
      solve();
   }

   /*
   * Compute stress from a polymer chain.
   */
   template <int D>
   void Polymer<D>::computeStress()
   {
      double prefactor;
      prefactor = 0;
     
      // Initialize stress_ to 0
      for (int i = 0; i < 6; ++i) {
        stress_[i] = 0.0;
      }

      for (int i = 0; i < nBlock(); ++i) {
         prefactor = exp(mu_)/length();
         block(i).computeStress(prefactor);
       
         for (int j=0; j < unitCellPtr_->nParameter() ; ++j){
            stress_[j] += block(i).stress(j);
         }
      }

   }


}
}
#endif
