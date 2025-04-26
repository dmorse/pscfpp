#ifndef RPC_BLOCK_TPP
#define RPC_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/cpu/FFT.h>
#include <pscf/chem/Edge.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <util/containers/DMatrix.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(nullptr),
      fftPtr_(nullptr),
      kMeshDimensions_(-1),
      kSize_(-1),
      ds_(-1.0),
      dsTarget_(-1.0),
      ns_(-1),
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

   /*
   * Store addresses of mesh, FFT and unit cell.
   */
   template <int D>
   void Block<D>::associate(Mesh<D> const & mesh,
                            FFT<D> const& fft,
                            UnitCell<D> const& cell)
   {
      // Preconditions
      UTIL_CHECK(!isAllocated_);

      // Set pointers to mesh and fft
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      unitCellPtr_ = &cell;

      hasExpKsq_ = false;
   }

   /*
   * Compute number of contour steps and allocate all memory.
   */
   template <int D>
   void Block<D>::allocate(double ds)
   {
      UTIL_CHECK(ds > 0.0);
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(fftPtr_);
      UTIL_CHECK(unitCellPtr_);
      UTIL_CHECK(mesh().size() > 1);
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(!isAllocated_);

      // Compute DFT grid dimensions kMeshDimensions_ and size kSize_
      FFT<D>::computeKMesh(mesh().dimensions(), kMeshDimensions_, kSize_);

      // Allocate work arrays for MDE solution
      expKsq_.allocate(kMeshDimensions_);
      expW_.allocate(mesh().dimensions());
      qr_.allocate(mesh().dimensions());
      qk_.allocate(mesh().dimensions());
      qr2_.allocate(mesh().dimensions());
      qk2_.allocate(mesh().dimensions());
      if (PolymerModel::isThread()) {
         expKsq2_.allocate(kMeshDimensions_);
         expW2_.allocate(mesh().dimensions());
      } else 
      if (PolymerModel::isBead()) {
         expWInv_.allocate(mesh().dimensions());
      }

      // Allocate work array for stress calculation
      dGsq_.allocate(kSize_, 6);

      // Allocate block concentration field
      cField().allocate(mesh().dimensions());

      dsTarget_ = ds;

      // Compute ns_ 
      if (PolymerModel::isThread()) {

         // Set contour length discretization for this block
         UTIL_CHECK(length() > 0.0);
         int tempNs;
         tempNs = floor( length()/(2.0 *ds) + 0.5 );
         if (tempNs == 0) {
            tempNs = 1;
         }
         ns_ = 2*tempNs + 1;
         ds_ = length()/double(ns_-1);

      } else
      if (PolymerModel::isBead()) {

         ds_ = 1.0;
         ns_ = nBead();
         if (!ownsVertex(0)) ++ns_;
         if (!ownsVertex(1)) ++ns_;

      }

      // Allocate memory for solutions to MDE (requires ns_)
      propagator(0).allocate(ns_, mesh());
      propagator(1).allocate(ns_, mesh());

      isAllocated_ = true;
      hasExpKsq_ = false;
   }

   /*
   * Set or reset the the block length.
   */
   template <int D>
   void Block<D>::setLength(double newLength)
   {
      UTIL_CHECK(PolymerModel::isThread());
      Edge::setLength(newLength);

      if (isAllocated_) {

         int oldNs = ns_;

         // Reset contour length discretization
         UTIL_CHECK(dsTarget_ > 0);
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
   * Mark data that depend on the unit cell parameters as invalid.
   */
   template <int D>
   void Block<D>::clearUnitCellData()
   {  hasExpKsq_ = false; }

   /*
   * Compute all elements of expKsq_ and expKsq2_ arrays
   */
   template <int D>
   void Block<D>::computeExpKsq()
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCellPtr_);
      UTIL_CHECK(unitCellPtr_->isInitialized());

      bool isThread = PolymerModel::isThread();
      double bSqFactor;
      if (isThread) {
         bSqFactor = -1.0*kuhn()*kuhn() * ds_ / 6.0;
      } else {
         bSqFactor = -1.0*kuhn()*kuhn() / 6.0;
      }


      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq, arg;
      int i;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank();
         G = iter.position();
         Gmin = shiftToMinimum(G, mesh().dimensions(), unitCell());
         Gsq = unitCell().ksq(Gmin);
         arg = Gsq*bSqFactor;
         expKsq_[i] = exp(arg);
         if (isThread) {
            expKsq2_[i] = exp(0.5*arg);
         }
      }

      hasExpKsq_ = true;
   }

   /*
   * Setup the the step algorithm for a specific field configuration.
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
      double arg;
      if (PolymerModel::isThread()) {
         double c = -0.5*ds_;
         for (int i = 0; i < nx; ++i) {
            arg = c*w[i];
            // UTIL_CHECK(std::abs(arg) < 0.5);
            expW_[i]  = exp(arg);
            expW2_[i] = exp(0.5*arg);
         }
      } else
      if (PolymerModel::isBead()) {
         for (int i = 0; i < nx; ++i) {
            arg = -w[i];
            // UTIL_CHECK(std::abs(arg) < 1.0);
            expW_[i]  = exp(arg);
            expWInv_[i] = 1.0/expW_[i];
         }
      }

      // Compute expKsq arrays if necessary
      if (!hasExpKsq_) {
         computeExpKsq();
      }

   }

   /*
   * Propagate solution by one step for the thread model.
   */
   template <int D>
   void Block<D>::stepThread(RField<D> const & q, RField<D>& qout)
   {
      UTIL_CHECK(PolymerModel::isThread());

      // Prereconditions on mesh and fft
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());

      // Internal preconditions
      UTIL_CHECK(isAllocated_);
      int nk = qk_.capacity();
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(expW2_.capacity() == nx);
      UTIL_CHECK(expKsq_.capacity() == nk);
      UTIL_CHECK(expKsq2_.capacity() == nk);
      UTIL_CHECK(hasExpKsq_);

      // Preconditions on parameters
      UTIL_CHECK(q.isAllocated());
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qout.isAllocated());
      UTIL_CHECK(qout.capacity() == nx);

      // Apply pseudo-spectral algorithm

      // Step by ds/2 for qr_, step by ds/4 for qr2_
      int i;
      for (i = 0; i < nx; ++i) {
         qr_[i] = q[i]*expW_[i];
         qr2_[i] = q[i]*expW2_[i];
      }
      fft().forwardTransform(qr_, qk_);
      fft().forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft().inverseTransformUnsafe(qk_, qr_); // overwrites qk_
      fft().inverseTransformUnsafe(qk2_, qr2_); // overwrites qk2_
      for (i = 0; i < nx; ++i) {
         qr_[i] = qr_[i]*expW_[i];
         qr2_[i] = qr2_[i]*expW_[i];
      }

      // Above, multiplying qr2_ by expW_ rather than expW2_ combines
      // required multiplications by expW2_ at the end of first half-step 
      // and at the beginning of the second.

      // Finish second half-step of ds/2 for qr2_
      fft().forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft().inverseTransformUnsafe(qk2_, qr2_); // overwrites qk2_
      for (i = 0; i < nx; ++i) {
         qr2_[i] = qr2_[i]*expW2_[i];
      }

      // Richardson extrapolation
      for (i = 0; i < nx; ++i) {
         qout[i] = (4.0*qr2_[i] - qr_[i])/3.0;
      }
   }

   /*
   * Apply one step of MDE solution for the bead model.
   */
   template <int D>
   void Block<D>::stepBead(RField<D> const & q, RField<D>& qout)
   {
      UTIL_CHECK(PolymerModel::isBead());
      stepBondBead(q, qout);
      stepFieldBead(qout);
   }

   /*
   * Apply the bond operator for the bead model.
   */
   template <int D>
   void Block<D>::stepBondBead(RField<D> const & q, RField<D>& qout)
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);

      // Prereconditions on mesh and fft
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qout.capacity() == nx);

      int nk = qk_.capacity();
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(expKsq_.capacity() == nk);

      // Apply bond operator
      fft().forwardTransform(q, qk_);
      for (int i = 0; i < nk; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
      }
      fft().inverseTransformUnsafe(qk_, qout); // destroys qk_

   }

   /*
   * Apply the local field operator for the bead model.
   */
   template <int D>
   void Block<D>::stepFieldBead(RField<D>& q)
   {
      // Preconditions 
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(q.capacity() == nx);

      // Apply field operator
      for (int i = 0; i < nx; ++i) {
         q[i] *= expW_[i];
      }
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentrationThread(double prefactor)
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

      // References to forward and reverse propagators
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
         RField<D> const & qf = p0.q(j);
         RField<D> const & qr = p1.q(ns_ - 1 - j);
         for (i = 0; i < nx; ++i) {
            //cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i] * 4.0;
            cField()[i] += qf[i] * qr[i] * 4.0;
         }
      }

      // Even indices
      for (j = 2; j < (ns_ -2); j += 2) {
         RField<D> const & qf = p0.q(j);
         RField<D> const & qr = p1.q(ns_ - 1 - j);
         for (i = 0; i < nx; ++i) {
            // cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i] * 2.0;
            cField()[i] += qf[i] * qr[i] * 2.0;
         }
      }

      // Normalize the integral
      prefactor *= ds_ / 3.0;
      for (i = 0; i < nx; ++i) {
         cField()[i] *= prefactor;
      }

   }

   /*
   * Calculate monomer concentration for this block, bead model.
   */
   template <int D>
   void Block<D>::computeConcentrationBead(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == nx);

      // Initialize cField to zero at all points
      int i;
      for (i = 0; i < nx; ++i) {
         cField()[i] = 0.0;
      }
   
      // References to forward and reverse propagators
      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      // Vertex 0 contribution (if owned by block)
      if (ownsVertex(0)) {
         RField<D> const & qf = p0.q(0);
         RField<D> const & qr = p1.q(ns_ - 1);
         for (i = 0; i < nx; ++i) {
            //cField()[i] += p0.q(0)[i]*p1.q(ns_ - 1)[i]*expWInv_[i];
            cField()[i] += qf[i] * qr[i] * expWInv_[i];
         }
      }

      // Internal beads (j = 1, ... , ns_ -2)
      int j;
      for (j = 1; j < (ns_ -1); ++j) {
         RField<D> const & qf = p0.q(j);
         RField<D> const & qr = p1.q(ns_ - 1 - j);
         for (i = 0; i < nx; ++i) {
            //cField()[i] += p0.q(j)[i] * p1.q(ns_ - 1 - j)[i]*expWInv_[i];
            cField()[i] += qf[i] * qr[i] * expWInv_[i];
         }
      }

      // Vertex 1 contribution (if owned by block)
      if (ownsVertex(1)) {
         RField<D> const & qf = p0.q(ns_ - 1);
         RField<D> const & qr = p1.q(0);
         for (i = 0; i < nx; ++i) {
            //cField()[i] += p0.q(ns_ -1)[i]*p1.q(0)[i]*expWInv_[i];
            cField()[i] += qf[i] * qr[i] * expWInv_[i];
         }
      }

      // Normalize the integral
      for (i = 0; i < nx; ++i) {
         cField()[i] *= prefactor;
      }

   }

   /*
   * Average of a product of complementary propagator slices.
   */
   template <int D>
   double 
   Block<D>::averageProduct(RField<D> const& q0, RField<D> const& q1)
   {
      int nx = mesh().size();
      UTIL_CHECK(nx == q0.capacity());
      UTIL_CHECK(nx == q1.capacity());

      double Q = 0.0; 
      for (int i = 0; i < nx; ++i) {
         Q += q0[i] * q1[i];
      }
      Q /= double(nx);
      return Q;
   }

   /*
   * Spatial integral of a product of complementary propagator slices.
   */
   template <int D>
   double 
   Block<D>::averageProductBead(RField<D> const& q0, RField<D> const& q1)
   {
      int nx = mesh().size();
      UTIL_CHECK(nx == q0.capacity());
      UTIL_CHECK(nx == q1.capacity());

      double Q = 0.0; 
      for (int i = 0; i < nx; ++i) {
         Q += q0[i] * q1[i] * expWInv_[i];
      }
      Q /= double(nx);
      return Q;
   }

   /*
   * Integrate to Stress exerted by the chain for this block
   */
   template <int D>
   void Block<D>::computeStressThread(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(PolymerModel::isThread());
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      computedGsq();
      stress_.clear();

      // Initialize work array and stress_ to zero at all points
      FSArray<double, 6> dQ;
      int nParam = unitCell().nParameter();
      for (int i = 0; i < nParam; ++i) {
         stress_.append(0.0);
         dQ.append(0.0);
      }

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      double dels, prod, increment;
      double bSq = kuhn()*kuhn()/6.0;
      int n, m;

      // Evaluate unnormalized integral over contour 
      for (int j = 0; j < ns_ ; ++j) {

         qr_ = p0.q(j);
         fft().forwardTransform(qr_, qk_);

         qr2_ = p1.q(ns_ - 1 - j);
         fft().forwardTransform(qr2_, qk2_);

         // Compute prefactor dels for Simpson's rule
         dels = ds_ / 3.0;
         if (j != 0 && j != ns_ - 1) {
            if (j % 2 == 0) {
               dels *= 2.0;
            } else {
               dels *= 4.0;
            }
         }

         // Loop over unit cell parameters
         for (n = 0; n < nParam ; ++n) {
            increment = 0.0;

            // Loop over wavevectors
            for (m = 0; m < kSize_ ; ++m) {
               prod = (qk2_[m][0] * qk_[m][0]) + (qk2_[m][1] * qk_[m][1]);
               prod *= dGsq_(m,n);
               increment += prod;
            }
            increment *= bSq * dels;
            dQ[n] = dQ[n] - increment;
         }

      }

      // Normalize
      for (int i = 0; i < nParam; ++i) {
         stress_[i] = stress_[i] - (dQ[i] * prefactor);
      }

   }

   /*
   * Compute contribution of this block to stress for bead model.
   */
   template <int D>
   void Block<D>::computeStressBead(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(PolymerModel::isBead());
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      computedGsq();
      stress_.clear();

      // Initialize dQ and stress_ to zero at all points
      FSArray<double, 6> dQ;
      int nParam = unitCell().nParameter();
      for (int i = 0; i < nParam; ++i) {
         dQ.append(0.0);
         stress_.append(0.0);
      }

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);
      double increment, prod;
      double bSq = kuhn()*kuhn()/6.0;

      // Loop over bonds in block
      for (int j = 0; j < ns_ - 1 ; ++j) {

         // Bead j, forward propagator
         qr_ = p0.q(j);
         fft().forwardTransform(qr_, qk_);

         // Bead j + 1, reverse propagator
         qr2_ = p1.q(ns_ - 2 - j);
         fft().forwardTransform(qr2_, qk2_);

         // Loop over unit cell parameters
         for (int n = 0; n < nParam ; ++n) {
            increment = 0.0;

            // Loop over wavevectors
            for (int m = 0; m < kSize_ ; ++m) {
               prod = (qk2_[m][0] * qk_[m][0]) + (qk2_[m][1] * qk_[m][1]);
               prod *= dGsq_(m, n)*expKsq_[m];
               increment += prod;
            }
            increment *= bSq;
            dQ[n] = dQ[n] - increment;
         }
      }

      // Normalize
      for (int i = 0; i < nParam; ++i) {
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

}
}
#endif
