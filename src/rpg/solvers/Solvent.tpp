#ifndef RPG_SOLVENT_TPP
#define RPG_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <prdc/cuda/resources.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rpg { 

   /*
   * Constructor.
   */
   template <int D>
   Solvent<D>::Solvent()
    : meshPtr_(nullptr)
   {  setClassName("Solvent"); }

   /*
   * Destructor.
   */
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

      // Evaluate unnormalized integral and Q
      double s = size();
      double Q = 0.0;

      // cField_ = exp(-size() * wField)
      VecOp::expVc(cField_, wField, -1.0*size());

      Q = Reduce::sum(cField_) / ((double) nx); // spatial average
      Q /= phiTot; // correct for partial occupation

      // Note: phiTot = 1.0 except in the case of a mask that confines
      // material to a fraction of the unit cell. 

      // Set q and compute mu or phi (depending on ensemble)
      Species::setQ(Q);

      // Normalize concentration 
      double prefactor = phi()/Q;
      VecOp::mulEqS(cField_, prefactor);
   }

}
}
#endif
