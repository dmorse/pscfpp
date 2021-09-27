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
      // Allocate concentration_ array
   }

   /*
   * Compute concentration, q, phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(WField const & wField)
   {
      // Setup s
   }

   /*
   * Compute solution to MDE and block concentrations.
   */ 
   template <int D>
   void Solvent<D>::compute(DArray<WField> const & wFields) 
   {  compute(wFields[monomerId()]); }

}
}
#endif
