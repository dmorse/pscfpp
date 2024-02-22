#ifndef RPG_SOLVENT_TPP
#define RPG_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rpg { 

   template <int D>
   Solvent<D>::Solvent()
   {  setClassName("Solvent"); }

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
      concField_.allocate(mesh.dimensions());
   }

   /*
   * Compute concentration, q, phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(RField<D> const & wField)
   {
      int nx = meshPtr_->size(); // Number of grid points

      // GPU Resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nx, nBlocks, nThreads);

      // Initialize concField_ to zero
      assignUniformReal<<<nBlocks, nThreads>>>(concField_.cField(), 0, nx);

      // Evaluate unnormalized integral and q_
      double s = size();
      q_ = 0.0;
      assignExp<<<nBlocks, nThreads>>>(concField_.cField(), wField.cField(), s, nx);
      q_ = (double)gpuSum(concField_.cField(),nx);
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
      scaleReal<<<nBlocks, nThreads>>>(concField_.cField(), prefactor, nx);
    
   }

}
}
#endif
