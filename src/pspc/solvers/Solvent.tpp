#ifndef PSPC_SOLVENT_TPP
#define PSPC_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Pspc { 

   /*
   * Constructor
   */
   template <int D>
   Solvent<D>::Solvent()
    : SolventDescriptor()
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
   void Solvent<D>::setDiscretization(Mesh<D> const & mesh)
   {
      meshPtr_ = &mesh;
      cField_.allocate(mesh.dimensions());
   }

   /*
   * Compute concentration, q, and phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(WField const & wField)
   {
      int nx = meshPtr_->size(); // Number of grid points

      // Initialize cField_ to zero
      for (int i = 0; i < nx; ++i) {
          cField_[i] = 0.0;
      }

      // Evaluate unnormalized integral and q_
      double s = size();
      q_ = 0.0;
      for (int i = 0; i < nx; ++i) {
          cField_[i] = exp(-s*wField[i]);
          q_ += cField_[i];
      }
      q_ = q_/double(nx);

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
      for (int i = 0; i < nx; ++i) {
          cField_[i] *= prefactor;
      }
    
   }

}
}
#endif
