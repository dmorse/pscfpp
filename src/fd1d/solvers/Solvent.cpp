#ifndef FD1D_SOLVENT_TPP
#define FD1D_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
      UTIL_CHECK(cField_.isAllocated());

      // Evaluate unnormalized concentration, Boltzmann weight
      int nx = domain().nx();
      double s = size();
      for (int i = 0; i < nx; ++i) {
          cField_[i] = exp(-s*wField[i]);
      }

      // Compute spatial average q_
      q_ = domain().spatialAverage(cField_);

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
