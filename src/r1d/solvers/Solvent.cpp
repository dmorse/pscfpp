/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <r1d/domain/Domain.h>

namespace Pscf {
namespace R1d { 

   /*
   * Constructor.
   */
   Solvent::Solvent() 
    : SolventSpecies()  
   {  setClassName("Solvent"); }

   /*
   * Destructor.
   */
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
   * Compute concentration, q, and phi or mu.
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

      // Compute spatial average q
      double Q = domain().spatialAverage(cField_);

      // Set q and compute mu or phi (depending on ensemble)
      Species::setQ(Q);

      #if 0
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
