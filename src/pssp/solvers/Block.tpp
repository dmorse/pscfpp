#ifndef PSSP_BLOCK_TPP
#define PSSP_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>

namespace Pscf { 
namespace Pssp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(0),
      kMeshDimensions_(0),
      ds_(0.0),
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
   void Block<D>::setDiscretization(double ds, Mesh<D>& mesh)
   {  
      UTIL_CHECK(mesh.size() > 1);
      UTIL_CHECK(ds > 0.0);

      // Set association to mesh
      meshPtr_ = &mesh;

      // Set contour length discretization
      ns_ = floor(length()/ds + 0.5) + 1;
      if (ns_%2 == 0) {
         ns_ += 1;
      }
      ds_ = length()/double(ns_ - 1);

      // Compute Fourier space kMeshDimensions_ 
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh.dimensions()[i];
         } else {
            kMeshDimensions_[i] = mesh.dimensions()[i]/2 + 1;
         }
      }

      // Allocate work arrays
      expKsq_.allocate(kMeshDimensions_);
      expW_.allocate(mesh.dimensions());
      qr_.allocate(mesh.dimensions());
      qk_.allocate(mesh.dimensions());

   }

   /*
   * Setup the contour length step algorithm.
   */
   template <int D>
   void 
   Block<D>::setupSolver(Block<D>::WField const& w, UnitCell<D>& unitCell)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      
      // Populate expW_
      int i;
      // std::cout << std::endl;
      for (i = 0; i < nx; ++i) {
         expW_[i] = exp(-0.5*w[i]*ds_);
         // std::cout << "i = " << i 
         //           << " expW_[i] = " << expW_[i]
         //          << std::endl;
      }

      MeshIterator<D> iter;
      IntVec<D> G;
      IntVec<D> Gmin;
      double Gsq;
      double factor = -1.0*kuhn()*kuhn()*ds_/6.0;
      // std::cout << "kDimensions = " << kMeshDimensions_ << std::endl;
      // std::cout << "factor      = " << factor << std::endl;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank(); 
         G = iter.position();
         Gmin = shiftToMinimum(G, mesh().dimensions(), unitCell);
         Gsq = unitCell.ksq(Gmin);
         expKsq_[i] = exp(Gsq*factor);
         //std::cout << i    << "  " 
         //          << Gmin << "  " 
         //          << Gsq  << "  "
         //          << expKsq_[i] << std::endl;
      }
      
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentration(double prefactor)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == nx) 

      // Initialize cField to zero at all points
      int i;
      for (i = 0; i < nx; ++i) {
         cField()[i] = 0.0;
      }

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

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
      #if 0
      #endif

   }

   /*
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::step(const QField& q, QField& qNew)
   {
      // Check real-space mesh sizes
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qNew.capacity() == nx);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);

      // Fourier-space mesh sizes
      int nk = qk_.capacity();
      UTIL_CHECK(expKsq_.capacity() == nk);

      // Apply pseudo-spectral algorithm
      int i;
      for (i = 0; i < nx; ++i) {
         qr_[i] = q[i]*expW_[i];
      }
      fft_.forwardTransform(qr_, qk_);
      for (i = 0; i < nk; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
      }
      fft_.inverseTransform(qk_, qr_);
      for (i = 0; i < nx; ++i) {
         qNew[i] = qr_[i]*expW_[i];
      }
   }

}
}
#endif
