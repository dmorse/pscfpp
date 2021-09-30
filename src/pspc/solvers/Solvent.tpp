#ifndef PSPC_SOLVENT_TPP
#define PSPC_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/mesh/Mesh.tpp>

namespace Pscf {
namespace Pspc { 

   template <int D>
   Solvent<D>::Solvent()
   {  setClassName("Solvent");}

   template <int D>
   Solvent<D>::~Solvent()
   {}

   template <int D>
   void Solvent<D>::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Species::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi;
   }

   template <int D>
   void Solvent<D>::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Species::Open);  
      mu_ = mu; 
   }

   template <int D>
   void Solvent<D>::setMesh(Mesh<D> const & mesh)
   {
      meshPtr_ = &mesh;
      cField_.allocate(mesh.dimensions());
   }

   /*
   * Compute concentration, q, phi or mu.
   */ 
   template <int D>
   void Solvent<D>::solve(WField const & wField)
   {
      int nx = mesh().size();
      int i;

      // Initialize to zero
      for (i = 0; i < nx; ++i) {
          cField_[i] = 0.0;
      }

      // Evaluate unnormalized integral and q_
      double s = size();
      q_ = 0.0;
      for (i = 0; i < nx; ++i) {
          cField_[i] = exp(-s*wField[i]);
          q_ += cField_[i];
      }
      q_ = q_/double(nx);

      // Compute mu_ or phi_ and prefactor
      double prefactor;
      if (ensemble_ == Species::Closed) {
         prefactor_ = phi_/q_;
         mu_ = log(prefactor);
      } else {
         prefactor_ = exp(mu_);
         phi_ = prefactor*q_;
      }

      // Normalize concentration 
      for (i = 0; i < nx; ++i) {
          cField_[i] *= prefactor;
      }
    
   }

}
}
#endif
