#ifndef PSPC_BLOCK_TPP
#define PSPC_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <util/containers/DMatrix.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(0),
      kMeshDimensions_(0),
      ds_(0.0),
      dsTarget_(0.0),
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
      dsTarget_ = ds;
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

      // Set number kSize of points in k-space mesh
      int kSize = 1;
      for (int i = 0; i < D; ++i) {
           kSize *= kMeshDimensions_[i];
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

      // Allocate work array for stress calculation
      dGsq_.allocate(kSize, 6);

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
   void Block<D>::setLength(double newLength)
   {
      BlockDescriptor::setLength(newLength);

      if (isAllocated_) { // if setDiscretization has already been called
         // Reset contour length discretization
         UTIL_CHECK(dsTarget_ > 0);
         int oldNs = ns_;
         int tempNs;
         tempNs = floor( length()/(2.0 *dsTarget_) + 0.5 );
         if (tempNs == 0) {
            tempNs = 1;
         }
         ns_ = 2*tempNs + 1;
         ds_ = length()/double(ns_-1);

         if (oldNs != ns_) {
            // If propagators are already allocated and ns_ has changed, 
            // reallocate memory for solutions to MDE
            propagator(0).reallocate(ns_);
            propagator(1).reallocate(ns_);
         }
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
      UTIL_CHECK(unitCellPtr_);
      UTIL_CHECK(unitCellPtr_->isInitialized());

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
   Block<D>::setupSolver(RField<D> const& w)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(isAllocated_);

      // Compute expW arrays
      for (int i = 0; i < nx; ++i) {

         // First, check that w[i]*ds_ is not unreasonably large:
         // (if this condition is not met, solution will have large
         // error, and user should consider using a smaller ds_)
         //UTIL_CHECK(std::abs(w[i]*ds_) < 1.0);

         // Calculate values
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
      UTIL_CHECK(cField().capacity() == nx);

      // Initialize cField to zero at all points
      int i;
      for (i = 0; i < nx; ++i) {
         cField()[i] = 0.0;
      }

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      // Evaluate unnormalized integral

      // Endpoint contributions
      for (i = 0; i < nx; ++i) {
         cField()[i] += p0.q(0)[i]*p1.q(ns_ - 1)[i];
         cField()[i] += p0.q(ns_ -1)[i]*p1.q(0)[i];
      }

      // Odd indices
      int j;
      for (j = 1; j < (ns_ -1); j += 2) {
         for (i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i] * 4.0;
         }
      }

      // Even indices
      for (j = 2; j < (ns_ -2); j += 2) {
         for (i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i] * 2.0;
         }
      }

      // Normalize the integral
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
      int nParam, c, m;

      normal = 3.0*6.0;

      int kSize_ = 1;
      for (int i = 0; i < D; ++i) {
           kSize_ *= kMeshDimensions_[i];
      }
      nParam = unitCell().nParameter();
      c = kSize_;

      FSArray<double, 6> dQ;

      // Initialize work array and stress_ to zero at all points
      int i;
      for (i = 0; i < nParam; ++i) {
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

         for (int n = 0; n < nParam ; ++n) {
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
      for (i = 0; i < nParam; ++i) {
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
   void Block<D>::step(RField<D> const & q, RField<D>& qNew)
   {
      // Internal prereconditions
      UTIL_CHECK(isAllocated_);
      int nx = mesh().size();
      int nk = qk_.capacity();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(expKsq_.capacity() == nk);
      UTIL_CHECK(hasExpKsq_);

      // Preconditions on parameters
      UTIL_CHECK(q.isAllocated());
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qNew.isAllocated());
      UTIL_CHECK(qNew.capacity() == nx);

      // Apply pseudo-spectral algorithm

      // Full step for ds, half-step for ds/2
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
         qr_[i] = qr_[i]*expW_[i];
         qr2_[i] = qr2_[i]*expW_[i];
      }

      // Finish second half-step for ds/2
      fft_.forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft_.inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         qr2_[i] = qr2_[i]*expW2_[i];
      }

      // Richardson extrapolation
      for (i = 0; i < nx; ++i) {
         qNew[i] = (4.0*qr2_[i] - qr_[i])/3.0;
      }
   }

}
}
#endif
