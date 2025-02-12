#ifndef RPG_SOLVENT_TPP
#define RPG_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <prdc/cuda/resources.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rpg { 

   template <int D>
   Solvent<D>::Solvent()
    : meshPtr_(nullptr)
   {  setClassName("Solvent"); }

   template <int D>
   Solvent<D>::~Solvent()
   {}

   /*
   * Compute concentration, q, phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(RField<D> const & wField, double phiTot)
   {
      int nx = meshPtr_->size(); // Number of grid points

      // Evaluate unnormalized integral and q_
      double s = size();
      q_ = 0.0;

      // cField_ = exp(-size() * wField)
      VecOp::expVc(cField_, wField, -1.0*size());

      q_ = Reduce::sum(cField_) / ((double) nx); // spatial average
      q_ /= phiTot; // correct for partial occupation

      // Note: phiTot = 1.0 except in the case of a mask that confines
      // material to a fraction of the unit cell. 

      // Compute mu_ or phi_ and prefactor
      double prefactor;
      if (ensemble_ == Species::Closed) {
         prefactor = phi_/q_;
         mu_ = log(prefactor);
      } else {
         prefactor = exp(mu_);
         phi_ = prefactor*q_;
      }

      // Normalize concentration 
      VecOp::mulEqS(cField_, prefactor);
   }

}
}
#endif
