#ifndef RPG_POLYMER_TPP
#define RPG_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pscf {
namespace Rpg { 

   template <int D>
   Polymer<D>::Polymer()
    : nParams_(0)
   {  setClassName("Polymer"); }

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
   * Store the number of lattice parameters in the unit cell.
   */ 
   template <int D>
   void Polymer<D>::setNParams(int nParams)
   {
      nParams_ = nParams;
   }

   /*
   * Compute solution to MDE and concentrations.
   */ 
   template <int D>
   void Polymer<D>::compute(DArray< RField<D> > const & wFields)
   {
      // Setup solvers for all blocks
      int monomerId;
      for (int j = 0; j < nBlock(); ++j) {
         monomerId = block(j).monomerId();
         block(j).setupSolver(wFields[monomerId]);
      }

      // Call generic solver() method base class template.
      // This solves the MDE for all propagators in a precalculated order
      solve();

      // Compute block concentration fields
      double prefactor;
      prefactor = phi() / ( q() * length() );
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }

   }

   /*
   * Compute stress from a polymer chain.
   */

   template <int D>
   void Polymer<D>::computeStress()
   {
      UTIL_CHECK(nParams_ > 0);
      
      double prefactor;
     
      // Initialize stress_ to 0
      for (int i = 0; i < nParams_; ++i) {
        stress_ [i] = 0.0;
      }

      for (int i = 0; i < nBlock(); ++i) {
         prefactor = exp(mu_)/length();
         block(i).computeStress(prefactor);
       
         for (int j=0; j < nParams_; ++j){
            stress_ [j] += block(i).stress(j);
            //std::cout<<"stress_[j] "<<stress_[j]<<std::endl;
         }
      }
   }

}
}
#endif
