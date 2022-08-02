#ifndef FD1D_SOLVENT_TPP
#define FD1D_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <fd1d/domain/Domain.h>

namespace Pscf {
namespace Fd1d { 

   Solvent::Solvent() 
    : SolventDescriptor()  
   {  setClassName("Solvent"); }

   Solvent::~Solvent()
   {}

   /*
   * Set association with Domain and allocate concentration field.
   */
   void Solvent::setDiscretization(Domain const & domain)
   {
      domainPtr_ = &domain;
      int nx = domain.nx();
      if (nx > 0) {
         cField_.allocate(nx);
      }
   }

   /*
   * Compute concentration, q, phi or mu.
   */ 
   void Solvent::compute(WField const & wField)
   {
      int nx = domainPtr_->nx();
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
         prefactor = phi_/q_;
         mu_ = log(prefactor);
      } else {
         prefactor = exp(mu_);
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
