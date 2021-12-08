#ifndef PSPC_BLOCK_TPP
#define PSPC_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/math/IntVec.h>
#include <util/containers/DMatrix.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(0),
      kMeshDimensions_(0),
      ds_(0.0),
      ns_(0),
      isAllocated_(false),
      hasExpKsq_(false)
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
   void Block<D>::setDiscretization(double ds, const Mesh<D>& mesh)
   {
      // Preconditions
      UTIL_CHECK(mesh.size() > 1);
      UTIL_CHECK(ds > 0.0);
      UTIL_CHECK(!isAllocated_);

      // Set contour length discretization
      int tempNs;
      tempNs = floor( length()/(2.0 *ds) + 0.5 );
      if (tempNs == 0) {
         tempNs = 1;
      }
      ns_ = 2*tempNs + 1;
      ds_ = length()/double(ns_-1);

      // Set association to mesh
      meshPtr_ = &mesh;

      fft_.setup(mesh.dimensions());

      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh.dimensions()[i];
         } else {
            kMeshDimensions_[i] = mesh.dimensions()[i]/2 + 1;
         }
      }

      // Set number kSize_ of points in k-space mesh
      int kSize_ = 1;
      for (int i = 0; i < D; ++i) {
           kSize_ *= kMeshDimensions_[i];
      }

      // Allocate work arrays for MDE solution
      expKsq_.allocate(kMeshDimensions_);
      expKsq2_.allocate(kMeshDimensions_);
      expW_.allocate(mesh.dimensions());
      expW2_.allocate(mesh.dimensions());
      qr_.allocate(mesh.dimensions());
      qk_.allocate(mesh.dimensions());
      qr2_.allocate(mesh.dimensions());
      qk2_.allocate(mesh.dimensions());
      qf_.allocate(mesh.dimensions());

      // Allocate work array for stress calculation
      dGsq_.allocate(kSize_, 6);

      // Allocate block concentration field
      cField().allocate(mesh.dimensions());

      // Allocate memory for solutions to MDE (requires ns_)
      propagator(0).allocate(ns_, mesh);
      propagator(1).allocate(ns_, mesh);
      
      isAllocated_ = true;
      hasExpKsq_ = false;
   }

   /*
   * Set or reset the the block length.
   */
   template <int D>
   void Block<D>::setLength(double length)
   {
      BlockDescriptor::setLength(length);
      if (isAllocated_) {
         UTIL_CHECK(ns_ > 1); 
         ds_ = length/double(ns_ - 1);
      }
      hasExpKsq_ = false;
   }

   /*
   * Set or reset the the block length.
   */
   template <int D>
   void Block<D>::setKuhn(double kuhn)
   {
      BlockTmpl< Propagator<D> >::setKuhn(kuhn);
      hasExpKsq_ = false;
   }

   /*
   * Setup data that depend on the unit cell parameters.
   */
   template <int D>
   void Block<D>::setupUnitCell(const UnitCell<D>& unitCell)
   {
      unitCellPtr_ = &unitCell;
      hasExpKsq_ = false;
   }

   template <int D>
   void Block<D>::computeExpKsq()
   {
      UTIL_CHECK(isAllocated_);

      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq;
      double factor = -1.0*kuhn()*kuhn()*ds_/6.0;
      int i;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank();
         G = iter.position();
         Gmin = shiftToMinimum(G, mesh().dimensions(), unitCell());
         Gsq = unitCell().ksq(Gmin);
         expKsq_[i] = exp(Gsq*factor);
         expKsq2_[i] = exp(Gsq*factor*0.5);
      }

      hasExpKsq_ = true;
   }

   /*
   * Setup the contour length step algorithm.
   */
   template <int D>
   void
   Block<D>::setupSolver(Block<D>::WField const& w)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(isAllocated_);

      // Compute expW arrays
      // double factor  = -0.5*ds_;
      // double factor2 = 0.5*factor;
      for (int i = 0; i < nx; ++i) {
         expW_[i] = exp(-0.5*w[i]*ds_);
         expW2_[i] = exp(-0.5*0.5*w[i]*ds_);
      }

      // Compute expKsq arrays if necessary
      if (!hasExpKsq_) {
         computeExpKsq();
      }

   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentration(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
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
      for(i = 0; i < nx; ++i) {
         cField()[i] += p0.q(0)[i]*p1.q(ns_ - 1)[i];
         cField()[i] += p0.q(ns_ -1)[i]*p1.q(0)[i];
      }

      //odd indices
      for(int j = 1; j < (ns_ -1); j += 2) {
         for(int i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i] * 4.0;
         }
      }

      //even indices
      for(int j = 2; j < (ns_ -2); j += 2) {
         for(int i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i] * 2.0;
         }
      }

      prefactor *= ds_ / 3.0;
      for (i = 0; i < nx; ++i) {
         cField()[i] *= prefactor;
      }

   }

   /*
   * Integrate to Stress exerted by the chain for this block
   */
   template <int D>
   void Block<D>::computeStress(double prefactor)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      stress_.clear();

      double dels, normal, increment;
      int r, c, m;

      normal = 3.0*6.0;

      int kSize_ = 1;
      for (int i = 0; i < D; ++i) {
           kSize_ *= kMeshDimensions_[i];
      }
      r = unitCell().nParameter();
      c = kSize_;

      FSArray<double, 6> dQ;

      // Initialize work array and stress_ to zero at all points
      int i;
      for (i = 0; i < r; ++i) {
         dQ.append(0.0);
         stress_.append(0.0);
      }

      computedGsq();

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      // Evaluate unnormalized integral
      for (int j = 0; j < ns_ ; ++j) {

         qr_ = p0.q(j);
         fft_.forwardTransform(qr_, qk_);

         qr2_ = p1.q(ns_ - 1 - j);
         fft_.forwardTransform(qr2_, qk2_);

         dels = ds_;

         if (j != 0 && j != ns_ - 1) {
            if (j % 2 == 0) {
               dels = dels*2.0;
            } else {
               dels = dels*4.0;
            }
         }

         for (int n = 0; n < r ; ++n) {
            increment = 0.0;

            for (m = 0; m < c ; ++m) {
               double prod = 0;
               prod = (qk2_[m][0] * qk_[m][0]) + (qk2_[m][1] * qk_[m][1]);
               prod *= dGsq_(m,n);
               increment += prod;
            }
            increment = (increment * kuhn() * kuhn() * dels)/normal;
            dQ[n] = dQ[n] - increment;
         }
      }

      // Normalize
      for (i = 0; i < r; ++i) {
         stress_[i] = stress_[i] - (dQ[i] * prefactor);
      }

   }

   /*
   * Compute dGsq_ array (derivatives of Gsq for all wavevectors)
   */
   template <int D>
   void Block<D>::computedGsq()
   {
      IntVec<D> temp;
      IntVec<D> vec;
      IntVec<D> Partner;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);

      for (int n = 0; n < unitCell().nParameter() ; ++n) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            temp = iter.position();
            vec = shiftToMinimum(temp, mesh().dimensions(), unitCell());
            dGsq_(iter.rank(), n) = unitCell().dksq(vec, n);
            for (int p = 0; p < D; ++p) {
               if (temp [p] != 0) {
                  Partner[p] = mesh().dimensions()[p] - temp[p];
               } else {
                  Partner[p] = 0;
               }
            }
            if (Partner[D-1] > kMeshDimensions_[D-1]) {
               dGsq_(iter.rank(), n) *= 2;
            }
         }
      }
   }

   /*
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::step(QField const & q, QField& qNew)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);

      // Check real space mesh sizes
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(q.isAllocated());
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qNew.isAllocated());
      UTIL_CHECK(qNew.capacity() == nx);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);

      // Check Fourier-space mesh sizes
      int nk = qk_.capacity();
      UTIL_CHECK(expKsq_.capacity() == nk);

      // Apply pseudo-spectral algorithm
      int i;
      for (i = 0; i < nx; ++i) {
         qr_[i] = q[i]*expW_[i];
         qr2_[i] = q[i]*expW2_[i];
      }
      fft_.forwardTransform(qr_, qk_);
      fft_.forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft_.inverseTransform(qk_, qr_);
      fft_.inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         qf_[i] = qr_[i]*expW_[i];
         qr2_[i] = qr2_[i]*expW_[i];
      }

      fft_.forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft_.inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         qr2_[i] = qr2_[i]*expW2_[i];
      }
      for (i = 0; i < nx; ++i) {
         qNew[i] = (4.0*qr2_[i] - qf_[i])/3.0;
      }
   }

}
}
#endif
