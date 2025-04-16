#ifndef RPC_SOLVENT_TPP
#define RPC_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rpc { 

   /*
   * Constructor
   */
   template <int D>
   Solvent<D>::Solvent()
    : SolventSpecies(),
      meshPtr_(nullptr)
   {  setClassName("Solvent"); }

   /*
   * Destructor
   */
   template <int D>
   Solvent<D>::~Solvent()
   {}

   /*
   * Create an association with a Mesh & allocate the concentration field.
   */
   template <int D>
   void Solvent<D>::associate(Mesh<D> const & mesh)
   {
      meshPtr_ = &mesh;
   }

   /*
   * Allocate the concentration field (cField).
   */
   template <int D>
   void Solvent<D>::allocate()
   {
      UTIL_CHECK(meshPtr_);
      cField_.allocate(meshPtr_->dimensions());
   }

   /*
   * Compute concentration, q, and phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(RField<D> const & wField, double phiTot)
   {
      int nx = meshPtr_->size(); // Number of grid points

      // Initialize cField_ to zero
      for (int i = 0; i < nx; ++i) {
          cField_[i] = 0.0;
      }

      // Evaluate unnormalized integral and Q
      double s = size();
      double Q = 0.0;
      for (int i = 0; i < nx; ++i) {
          cField_[i] = exp(-s*wField[i]);
          Q += cField_[i];
      }
      Q = Q/double(nx);  // spatial average
      Q = Q/phiTot;      // correct for partial occupation

      // Note: phiTot = 1.0 except in the case of a mask that confines
      // material to a fraction of the unit cell. 

      // Set q and compute mu or phi 
      Species::setQ(Q);

      #if 0
      q_ = Q;
      double prefactor;
      if (ensemble() == Species::Closed) {
         prefactor = phi_/q_;
         mu_ = log(prefactor);
      } else {
         prefactor = exp(mu_);
         phi_ = prefactor*q_;
      }
      #endif

      // Normalize concentration 
      double prefactor = phi()/Q;
      for (int i = 0; i < nx; ++i) {
          cField_[i] *= prefactor;
      }
 
   }

}
}
#endif
