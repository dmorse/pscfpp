#ifndef CPC_BLOCK_TPP
#define CPC_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include "Propagator.h"

#include <prdc/field/cpu/WaveList.h>
#include <prdc/field/cpu/FFT.h>
#include <pscf/cpu/complex.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/chem/Edge.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/solvers/BlockTmpl.tpp>

#include <util/containers/DMatrix.h>
#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>

#include <cmath>

namespace Pscf {
namespace Cpc {

   // Namespaces from which names may be used without qualification
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
      unitCellPtr_(nullptr),
      waveListPtr_(nullptr),
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
                            UnitCell<D> const& cell,
                            WaveList<D>& wavelist)
   {
      // Preconditions
      UTIL_CHECK(!isAllocated_);

      // Set pointers to mesh and fft
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      unitCellPtr_ = &cell;
      waveListPtr_ = &wavelist;

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
      UTIL_CHECK(waveListPtr_);
      UTIL_CHECK(mesh().size() > 1);
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(!isAllocated_);

      // Allocate work arrays for MDE solution
      expKsq_.allocate(mesh().dimensions());
      expKsq2_.allocate(mesh().dimensions());
      expW_.allocate(mesh().dimensions());
      qr_.allocate(mesh().dimensions());
      qr2_.allocate(mesh().dimensions());
      qk_.allocate(mesh().dimensions());
      qk2_.allocate(mesh().dimensions());
      if (PolymerModel::isThread()) {
         expW2_.allocate(mesh().dimensions());
      } else {
         expWInv_.allocate(mesh().dimensions());
      }

      // Allocate block concentration field
      cField().allocate(mesh().dimensions());

      dsTarget_ = ds;

      // Compute ns_
      if (PolymerModel::isThread()) {
         // Set contour length discretization for this block
         UTIL_CHECK(length() > 0.0);
         int tempNs;
         tempNs = floor(length() / (2.0 *ds) + 0.5);
         if (tempNs == 0) {
            tempNs = 1;                     // ensure ns_ >= 3
         }
         ns_ = 2*tempNs + 1;
         ds_ = length()/double(ns_-1);
      } else
      if (PolymerModel::isBead()) {
         ds_ = 1.0;
         ns_ = nBead() + 2;
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

         // Reallocate propagators if ns_ has changed
         if (oldNs != ns_) {
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
      Base::setKuhn(kuhn);
      hasExpKsq_ = false;
   }

   /*
   * Clear all internal data that depends on the unit cell parameters.
   */
   template <int D>
   void Block<D>::clearUnitCellData()
   {
      hasExpKsq_ = false;
   }

   /*
   * Compute all elements of expKsq_ and expKsq2_ arrays
   */
   template <int D>
   void Block<D>::computeExpKsq()
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(waveListPtr_);
      UTIL_CHECK(waveListPtr_->isAllocated());

      // Calculate KSq if necessary
      if (!waveListPtr_->hasKSq()) {
         waveListPtr_->computeKSq();
      }
      RField<D> const & kSq = waveListPtr_->kSq();
      UTIL_CHECK(kSq.capacity() == mesh().size());

      // Compute bSqFactor = b*b*ds/6
      bool isThread = PolymerModel::isThread();
      double bSqFactor = -1.0*kuhn()*kuhn() / 6.0;
      if (isThread) {
         bSqFactor *= ds_;
      }

      // Compute expKsq arrays
      MeshIterator<D> iter;
      iter.setDimensions(mesh().dimensions());
      double arg;
      int i;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank();
         arg = kSq[i]*bSqFactor;
         expKsq_[i] = exp(arg);
         expKsq2_[i] = exp(0.5*arg);
      }

      hasExpKsq_ = true;
   }

   /*
   * Setup the the step algorithm for a specific field configuration.
   */
   template <int D>
   void
   Block<D>::setupSolver(CField<D> const& w)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(isAllocated_);

      // Compute expW arrays
      fftw_complex arg, arg2;
      if (PolymerModel::isThread()) {
         double c = -0.5*ds_;
         double c2 = -0.25*ds_;
         for (int i = 0; i < nx; ++i) {

            mul(arg, w[i], c);
            mul(arg2, w[i], c2);
            assignExp(expW_[i], arg);
            assignExp(expW2_[i], arg2);

            //arg = c*w[i];
            //expW_[i]  = exp(arg);
            //expW2_[i] = exp(0.5*arg);

         }
      } else
      if (PolymerModel::isBead()) {
         for (int i = 0; i < nx; ++i) {

            mul(arg, w[i], -1.0);
            assignExp(expW_[i], arg);
            inverse(expWInv_[i], expW_[i]);

            //arg = -w[i];
            //expW_[i]  = exp(arg);
            //expWInv_[i] = 1.0/expW_[i];

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
   void Block<D>::stepThread(CField<D> const & q, CField<D>& qout) const
   {
      UTIL_CHECK(PolymerModel::isThread());

      // Preconditions on mesh and fft
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());

      // Internal preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(qk_.capacity() == nx);
      UTIL_CHECK(expKsq_.capacity() == nx);
      UTIL_CHECK(expKsq2_.capacity() == nx);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(expW2_.capacity() == nx);
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

         mul(qr_[i],  q[i], expW_[i] );
         mul(qr2_[i], q[i], expW2_[i]);

         // qr_[i] = q[i]*expW_[i];
         // qr2_[i] = q[i]*expW2_[i];
      }
      fft().forwardTransform(qr_, qk_);
      fft().forwardTransform(qr2_, qk2_);
      for (i = 0; i < nx; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft().inverseTransform(qk_, qr_);
      fft().inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         mulEq(qr_[i],  expW_[i]);
         mulEq(qr2_[i], expW_[i]);
         // qr_[i] = qr_[i]*expW_[i];
         // qr2_[i] = qr2_[i]*expW_[i];
      }

      // Note: Above, multiplying qr2_ by expW_, rather than by expW2_,
      // combines required multiplications by expW2_ at the end of first
      // half-step and at the beginning of the second.

      // Finish second half-step of ds/2 for qr2_
      fft().forwardTransform(qr2_, qk2_);
      for (i = 0; i < nx; ++i) {
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft().inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         mulEq(qr2_[i], expW2_[i]);
         //qr2_[i] = qr2_[i]*expW2_[i];
      }

      // Richardson extrapolation
      for (i = 0; i < nx; ++i) {
         qout[i][0] = (4.0*qr2_[i][0] - qr_[i][0])/3.0;
         qout[i][1] = (4.0*qr2_[i][1] - qr_[i][1])/3.0;
         //qout[i] = (4.0*qr2_[i] - qr_[i])/3.0;
      }

   }

   /*
   * Apply one step of MDE solution for the bead model.
   */
   template <int D>
   void Block<D>::stepBead(CField<D> const & q, CField<D>& qout) const
   {
      UTIL_CHECK(PolymerModel::isBead());
      stepBondBead(q, qout);
      stepFieldBead(qout);
   }

   /*
   * Apply the bond operator for the bead model.
   */
   template <int D>
   void Block<D>::stepBondBead(CField<D> const & q,
                               CField<D> & qout) const
   {
      // Prereconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qout.capacity() == nx);
      UTIL_CHECK(expKsq_.capacity() == nx);

      // Apply bond operator
      fft().forwardTransform(q, qk_);
      for (int i = 0; i < nx; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
      }
      fft().inverseTransform(qk_, qout);
   }

   /*
   * Apply the half-bond operator for the bead model.
   */
   template <int D>
   void Block<D>::stepHalfBondBead(CField<D> const & q,
                                   CField<D> & qout) const
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qout.capacity() == nx);
      UTIL_CHECK(qk_.capacity() == nx);
      UTIL_CHECK(expKsq2_.capacity() == nx);

      // Apply bond operator
      fft().forwardTransform(q, qk_);
      for (int i = 0; i < nx; ++i) {
         qk_[i][0] *= expKsq2_[i];
         qk_[i][1] *= expKsq2_[i];
      }
      fft().inverseTransform(qk_, qout); // destroys qk_
   }

   /*
   * Apply the local field operator for the bead model.
   */
   template <int D>
   void Block<D>::stepFieldBead(CField<D>& q) const
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(q.capacity() == nx);

      // Apply field operator
      for (int i = 0; i < nx; ++i) {
         mulEq(q[i], expW_[i]);
         //q[i] *= expW_[i];
      }
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void
   Block<D>::computeConcentrationThread(fftw_complex const & prefactor)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isSolved());
      UTIL_CHECK(propagator(1).isSolved());
      UTIL_CHECK(cField().capacity() == nx);

      CField<D> & c = cField();
      double d;
      int i, j;

      // Initialize cField to zero at all points
      d = 0.0;
      for (i = 0; i < nx; ++i) {
         assign(c[i], d);
      }

      // References to forward and reverse propagators
      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);
      fftw_complex z;

      // Evaluate unnormalized integral

      // Initial (head) endpoint contribution
      {
         CField<D> const & hf = p0.q(0);
         CField<D> const & hr = p1.q(ns_ - 1);
         for (i = 0; i < nx; ++i) {
            mul(z, hf[i], hr[i]);
            addEq(c[i], z);
         }
      }

      // Final (tail) endpoint contribution
      {
         CField<D> const & tf = p0.q(ns_ - 1);
         CField<D> const & tr = p1.q(0);
         for (i = 0; i < nx; ++i) {
            mul(z, tf[i], tr[i]);
            addEq(c[i], z);
         }
      }

      // Odd indices
      d = 4.0;
      for (j = 1; j < (ns_ - 1); j += 2) {
         CField<D> const & qf = p0.q(j);
         CField<D> const & qr = p1.q(ns_ - 1 - j);
         for (i = 0; i < nx; ++i) {
            mul(z, qf[i], qr[i]);
            mulEq(z, d);
            addEq(c[i], z);
            //c[i] += qf[i] * qr[i] * 4.0;
         }
      }

      // Even indices
      d = 2.0;
      for (j = 2; j < (ns_ - 2); j += 2) {
         CField<D> const & qf = p0.q(j);
         CField<D> const & qr = p1.q(ns_ - 1 - j);
         for (i = 0; i < nx; ++i) {
            mul(z, qf[i], qr[i]);
            mulEq(z, d);
            addEq(c[i], z);
            // c[i] += qf[i] * qr[i] * 2.0;
         }
      }

      // Normalize the integral
      d = ds_ / 3.0;
      mul(z, prefactor, d);
      // z = prefactor * ds_ / 3.0
      for (i = 0; i < nx; ++i) {
         mulEq(c[i], z);
         // c[i] *= z;
      }

   }

   /*
   * Calculate monomer concentration for this block, bead model.
   */
   template <int D>
   void Block<D>::computeConcentrationBead(fftw_complex const & prefactor)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isSolved());
      UTIL_CHECK(propagator(1).isSolved());
      UTIL_CHECK(cField().capacity() == nx);

      CField<D> & c = cField();
      fftw_complex z;
      double d;
      int i, j;

      // Initialize cField to zero at all points
      d = 0.0;
      for (i = 0; i < nx; ++i) {
         assign(c[i], d);
         // c[i] = 0.0;
      }

      // References to forward and reverse propagators
      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      // Sum over interior beads (j = 1, ... , ns_ -2)
      for (j = 1; j < (ns_ -1); ++j) {
         CField<D> const & qf = p0.q(j);
         CField<D> const & qr = p1.q(ns_ - 1 - j);
         for (i = 0; i < nx; ++i) {
            mul(z, qf[i], qr[i]);
            mulEq(z, expWInv_[i]);
            addEq(c[i], z);
            //c[i] += qf[i] * qr[i] * expWInv_[i];
         }
      }

      // Note: Slices j = 0 and j = ns_ - 1 are phantom vertices

      // Normalize the integral
      for (i = 0; i < nx; ++i) {
         mulEq(c[i], prefactor);
         //c[i] *= prefactor;
      }
   }

}
}
#endif
