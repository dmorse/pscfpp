#ifndef PSSP_BLOCK_TPP
#define PSSP_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf { 
namespace Pssp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : ds_(0.0),
      ns_(0)
   {
      propagator(0).setBlock(*this);
      propagator(1).setBlock(*this);
   }

   /*
   * Destructor.
   */
   template <int D>
   Block<D>::~Block()
   {}

   template <int D>
   //void Block<D>::setDiscretization(Domain const & domain, double ds)
   void Block<D>::setDiscretization(double ds, Mesh<D>& mesh)
   {  
      UTIL_CHECK(length() > 0);
      //UTIL_CHECK(domain.nx() > 1);
      UTIL_CHECK(ds > 0.0);

      // Set association to spatial domain
      // domainPtr_ = &domain;

      // Set contour length discretization
      ns_ = floor(length()/ds + 0.5) + 1;
      if (ns_%2 == 0) {
         ns_ += 1;
      }
      ds_ = length()/double(ns_ - 1);

      // Allocate all required memory
   }

   /*
   * Setup the contour length step algorithm.
   */
   template <int D>
   void Block<D>::setupSolver(Block<D>::WField const& w, UnitCell<D>& unitCell)
   {
      // Preconditions
      // UTIL_CHECK(domainPtr_);
      // int nx = domain().nx();

   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentration(double prefactor)
   {
      // Preconditions
      // UTIL_CHECK(domain().nx() > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      // UTIL_CHECK(cField().capacity() == domain().nx()) 

      #if 0
      // Initialize cField to zero at all points
      int i;
      //int nx = domain().nx();
      for (i = 0; i < nx; ++i) {
         cField()[i] = 0.0;
      }

      Propagator const & p0 = propagator(0);
      Propagator const & p1 = propagator(1);

      // Evaluate unnormalized integral
      for (i = 0; i < nx; ++i) {
         cField()[i] += 0.5*p0.q(0)[i]*p1.q(ns_ - 1)[i];
      }
      for (int j = 1; j < ns_ - 1; ++j) {
         for (i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i]*p1.q(ns_ - 1 - j)[i];
         }
      }
      for (i = 0; i < nx; ++i) {
         cField()[i] += 0.5*p0.q(ns_ - 1)[i]*p1.q(0)[i];
      }

      // Normalize
      prefactor *= ds_;
      for (i = 0; i < nx; ++i) {
         cField()[i] *= prefactor;
      }
      #endif

   }

   /*
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::step(const QField& q, QField& qNew)
   {
      //int nx = domain().nx();
   }

}
}
#endif
