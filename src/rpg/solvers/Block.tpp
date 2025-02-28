#ifndef RPG_BLOCK_TPP
#define RPG_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <rpg/solvers/WaveList.h>

#include <prdc/cuda/resources.h>
#include <prdc/cuda/FFT.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   // CUDA kernels:
   // (defined in anonymous namespace, used only in this file)

   namespace {

      /*
      * Element-wise calculation of a = real(b * conj(c) * d), CUDA kernel
      */
      __global__ void _realMulVConjVV(cudaReal* a, cudaComplex const * b,
                                      cudaComplex const * c,
                                      cudaReal const * d, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaComplex bt, ct;
         for (int i = startID; i < n; i += nThreads) {
            // Load complex numbers from global memory to local
            // (this way the accessed memory is contiguous, rather than
            // accessing the x or y elements individually, which are not
            // contiguous)
            bt = b[i];
            ct = c[i];

            // Perform calculation
            a[i] = ((bt.x * ct.x) + (bt.y * ct.y)) * d[i];
         }
      }

      /*
      * Performs qNew = (4*(qr2 * expW2)-qr)/3 elementwise, CUDA kernel
      */
      __global__ void _richardsonEx(cudaReal* qNew, 
                                    cudaReal const * qr,
                                    cudaReal const * qr2,
                                    cudaReal const * expW2, const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         cudaReal q2;
         for (int i = startID; i < n; i += nThreads) {
            q2 = qr2[i] * expW2[i];
            qNew[i] = (4.0 * q2 - qr[i]) / 3.0;
         }
      }

      /*
      * Performs a[i] += b[i] * c[i] * d. CUDA kernel
      */
      __global__ void _addEqMulVVc(cudaReal* a, 
                                   cudaReal const * b,
                                   cudaReal const * c, 
                                   const cudaReal d,
                                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b[i] * c[i] * d;
         }
      }

      /*
      * Element-wise calculation of a[i] += b[i]*c[i]* d[i], CUDA kernel
      */
      __global__ void _addEqMulVVV(cudaReal* a, 
                                   cudaReal const * b,
                                   cudaReal const * c,
                                   cudaReal const * d, 
                                   const int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] += b[i]*c[i]*d[i];
         }
      }

   }

   // CUDA kernel wrappers:

   /*
   * Element-wise calculation of a = real(b * conj(c) * d), kernel wrapper
   * 
   * \param a  output array (real)
   * \param b  input array 1 (complex)
   * \param c  input array 2 (complex)
   * \param d  input array 3 (real)
   */
   void realMulVConjVV(DeviceArray<cudaReal>& a,
                       DeviceArray<cudaComplex> const & b,
                       DeviceArray<cudaComplex> const & c,
                       DeviceArray<cudaReal> const & d)
   {
      int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      UTIL_CHECK(d.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _realMulVConjVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(),
                                             c.cArray(), d.cArray(), n);
   }

   /*
   * Performs qNew = (4 * (qr2 * expW2) - qr) / 3 elementwise, kernel wrapper
   * 
   * \param qNew  output array (a propagator slice)
   * \param qr  input array 1 (a propagator slice)
   * \param qr2  input array 2 (a propagator slice)
   * \param expW2  input array 3 (exp(-W[i]*ds/4) array)
   */
   void richardsonEx(DeviceArray<cudaReal>& qNew,
                     DeviceArray<cudaReal> const & qr,
                     DeviceArray<cudaReal> const & qr2,
                     DeviceArray<cudaReal> const & expW2)
   {
      int n = qNew.capacity();
      UTIL_CHECK(qr.capacity() == n);
      UTIL_CHECK(qr2.capacity() == n);
      UTIL_CHECK(expW2.capacity() == n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _richardsonEx<<<nBlocks, nThreads>>>(qNew.cArray(), qr.cArray(),
                                           qr2.cArray(), expW2.cArray(), n);
   }

   /*
   * Performs a[i] += b[i] * c[i] * d, kernel wrapper
   * 
   * \param a  output array
   * \param b  input array 1
   * \param c  input array 2
   * \param d  input scalar
   */
   void addEqMulVVc(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
                    DeviceArray<cudaReal> const & c, cudaReal const d)
   {
      int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqMulVVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(),
                                          c.cArray(), d, n);
   }

   /*
   * Element-wise calculation of a[i] = b[i]*c[i]*d[i], kernel wrapper
   */
   void addEqMulVVV(DeviceArray<cudaReal>& a,
                    DeviceArray<cudaReal> const & b,
                    DeviceArray<cudaReal> const & c,
                    DeviceArray<cudaReal> const & d)
   {
      int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      UTIL_CHECK(d.capacity() >= n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqMulVVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(),
                                          c.cArray(), d.cArray(), n);
   }

   // Block<D> member functions:

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(nullptr),
      fftPtr_(nullptr),
      unitCellPtr_(nullptr),
      waveListPtr_(nullptr),
      kMeshDimensions_(0),
      kSize_(0),
      ds_(0.0),
      dsTarget_(0.0),
      ns_(0),
      isAllocated_(false),
      hasExpKsq_(false),
      useBatchedFFT_(true),
      nParams_(0)
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
   void Block<D>::associate(Mesh<D> const & mesh, FFT<D> const & fft,
                            UnitCell<D> const & cell, WaveList<D>& wavelist)
   {
      UTIL_CHECK(!isAllocated_);
      UTIL_CHECK(mesh.size() > 1);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(mesh.dimensions() == fft.meshDimensions());

      nParams_ = cell.nParameter();
      UTIL_CHECK(nParams_ > 0);

      // Store pointers to associated objects
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      unitCellPtr_ = &cell;
      waveListPtr_ = &wavelist;

      hasExpKsq_ = false;
   }

   template <int D>
   void Block<D>::allocate(double ds, bool useBatchedFFT)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(unitCellPtr_);
      UTIL_CHECK(ds > 0.0);
      UTIL_CHECK(!isAllocated_);

      // Store useBatchedFFT
      useBatchedFFT_ = useBatchedFFT;

      // Compute Fourier space kMeshDimensions_ and kSize_
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh().dimensions()[i];
         } else {
            kMeshDimensions_[i] = mesh().dimensions()[i]/2 + 1;
         }
         kSize_ *= kMeshDimensions_[i];
      }

      // Allocate work arrays
      expW_.allocate(mesh().dimensions());
      expKsq_.allocate(kMeshDimensions_);
      if (PolymerModel::isThread()) {
         expW2_.allocate(mesh().dimensions());
         expKsq2_.allocate(kMeshDimensions_);
         qrPair_.allocate(2 * mesh().size());
         qkPair_.allocate(2 * kSize_);
      } else 
      if (PolymerModel::isBead()) {
         expWInv_.allocate(mesh().dimensions());
         qk_.allocate(mesh().dimensions());
      }

      // Allocate space for block monomer concentration
      cField().allocate(mesh().dimensions());

      // Compute ns_
      dsTarget_ = ds;
      if (PolymerModel::isThread()) {

         // Set contour length discretization for this block
         UTIL_CHECK(length() > 0.0);
         int tempNs;
         tempNs = floor(length() / (2.0 * ds) + 0.5);
         if (tempNs == 0) {
            tempNs = 1; // ensure at least 3 contour steps per chain
         }
         ns_ = 2*tempNs + 1;
         ds_ = length()/double(ns_ - 1);

      } else
      if (PolymerModel::isBead()) {

         ds_ = ds;
         ns_ = nBead();
         if (!ownsVertex(0)) ++ns_;
         if (!ownsVertex(1)) ++ns_;

      }

      // Allocate memory for solutions of MDE (requires ns_)
      propagator(0).allocate(ns_, mesh());
      propagator(1).allocate(ns_, mesh());

      // Setup fftBatchedPair_
      UTIL_CHECK(!fftBatchedPair_.isSetup());
      fftBatchedPair_.setup(mesh().dimensions(), 2);

      // Setup batched data used for stress calculation
      if (useBatchedFFT_) {
         UTIL_CHECK(!fftBatchedAll_.isSetup());
         fftBatchedAll_.setup(mesh().dimensions(), ns_);
         q0kBatched_.allocate(ns_ * kSize_);
         q1kBatched_.allocate(ns_ * kSize_);
      }

      isAllocated_ = true;
      hasExpKsq_ = false;
   }

   /*
   * Clear all internal data that depends on lattice parameters.
   */
   template <int D>
   void Block<D>::clearUnitCellData()
   {
      UTIL_CHECK(unitCellPtr_);
      UTIL_CHECK(nParams_ == unitCell().nParameter());
      hasExpKsq_ = false;
   }

   /*
   * Set or reset the the block length.
   */
   template <int D>
   void Block<D>::setLength(double newLength)
   {
      // Precondition
      UTIL_CHECK(PolymerModel::isThread());

      BlockDescriptor::setLength(newLength);

      if (isAllocated_) {
         // Reset contour length discretization
         UTIL_CHECK(dsTarget_ > 0);
         int oldNs = ns_;
         int tempNs;
         tempNs = floor(length() / (2.0 * dsTarget_) + 0.5);
         if (tempNs == 0) {
            tempNs = 1; // ensure at least 3 contour steps per chain
         }
         ns_ = 2*tempNs + 1;
         ds_ = length()/double(ns_-1);

         if (oldNs != ns_) {
            // If propagators are already allocated and ns_ has changed,
            // reallocate memory for solutions to MDE
            propagator(0).reallocate(ns_);
            propagator(1).reallocate(ns_);

            // If using batched FFTs, resize arrays and change batch size
            if (useBatchedFFT_) {
               UTIL_CHECK(fftBatchedAll_.isSetup());
               q0kBatched_.deallocate();
               q1kBatched_.deallocate();
               q0kBatched_.allocate(ns_ * kSize_);
               q1kBatched_.allocate(ns_ * kSize_);
               fftBatchedAll_.resetBatchSize(ns_);
            }
         }
      }

      hasExpKsq_ = false;
   }

   /*
   * Set or reset monomer statistical segment length.
   */
   template <int D>
   void Block<D>::setKuhn(double kuhn)
   {
      BlockTmpl< Propagator<D> >::setKuhn(kuhn);
      hasExpKsq_ = false;
   }

   template <int D>
   void Block<D>::computeExpKsq()
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCellPtr_);
      UTIL_CHECK(unitCellPtr_->isInitialized());

      // Calculate kSq if necessary
      if (!waveListPtr_->hasKSq()) {
         waveListPtr_->computeKSq();
      }

      // Calculate expKsq values on device
      if (PolymerModel::isThread()) {
         double bSqFactor = -1.0 * kuhn() * kuhn() * ds_ / 6.0;
         VecOp::expVc(expKsq_, waveListPtr_->kSq(), bSqFactor);
         VecOp::expVc(expKsq2_, waveListPtr_->kSq(), bSqFactor / 2.0);
      } else {
         double bSqFactor = -1.0 * kuhn() * kuhn() / 6.0;
         VecOp::expVc(expKsq_, waveListPtr_->kSq(), bSqFactor);
      }

      hasExpKsq_ = true;
   }

   /*
   * Setup the contour length step algorithm.
   */
   template <int D>
   void
   Block<D>::setupSolver(RField<D> const & w)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(isAllocated_);

      // Populate expW_
      if (PolymerModel::isThread()) {
         VecOp::expVc(expW_, w, -0.5 * ds_);
         VecOp::expVc(expW2_, w, -0.25 * ds_);
      } else {
         VecOp::expVc(expW_, w, -1.0);
         VecOp::divSV(expWInv_, 1.0, expW_);
      }

      // Compute expKsq arrays if necessary
      if (!hasExpKsq_) {
         computeExpKsq();
      }

   }

   /*
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::stepThread(RField<D> const & qin, RField<D>& qout)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);

      // Check real-space mesh sizes
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(fft().meshDimensions() == mesh().dimensions());
      UTIL_CHECK(qrPair_.capacity() == nx * 2);
      UTIL_CHECK(qkPair_.capacity() == kSize_ * 2);
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(expKsq_.capacity() == kSize_);
      UTIL_CHECK(fftBatchedPair_.isSetup());

      // Set up associated workspace fields slices
      RField<D> qr, qr2;
      RFieldDft<D> qk, qk2;
      qr.associate(qrPair_, 0, mesh().dimensions());
      qr2.associate(qrPair_, nx, mesh().dimensions());
      qk.associate(qkPair_, 0, mesh().dimensions());
      qk2.associate(qkPair_, kSize_, mesh().dimensions());

      // Apply pseudo-spectral algorithm
      VecOp::mulVVPair(qr, qr2, expW_, expW2_, qin); // qr = expW*q, qr2 = expW2*q
      fftBatchedPair_.forwardTransform(qrPair_, qkPair_); // to Fourier space
      VecOp::mulEqV(qk, expKsq_); // qk *= expKsq
      VecOp::mulEqV(qk2, expKsq2_); // qk2 *= expKsq2
      fftBatchedPair_.inverseTransformUnsafe(qkPair_, qrPair_); // to real space
      VecOp::mulEqVPair(qr, qr2, expW_); // qr *= expW, qr2 *= expW
      fft().forwardTransform(qr2, qk2); // to Fourier space, only qr2
      VecOp::mulEqV(qk2, expKsq2_); // qk2 *= expKsq2
      fft().inverseTransformUnsafe(qk2, qr2); // to real space, only qr2
      richardsonEx(qout, qr, qr2, expW2_); // qout=(4*(qr2*expW2)-qr)/3
   }

   /*
   * Apply one step of the MDE solution for the bead model. 
   */
   template <int D>
   void Block<D>::stepBead(RField<D> const & qin, RField<D>& qout)
   {
      stepBondBead(qin, qout);
      stepFieldBead(qout);
   }

   /*
   * Apply the bond operator for the bead model.
   */
   template <int D>
   void Block<D>::stepBondBead(RField<D> const & qin, RField<D>& qout)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);
      UTIL_CHECK(fft().isSetup());
      UTIL_CHECK(fft().meshDimensions() == mesh().dimensions());
      UTIL_CHECK(qk_.isAllocated());

      // Check mesh sizes
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(qin.capacity() == nx);
      UTIL_CHECK(qout.capacity() == nx);
      UTIL_CHECK(qk_.capacity() == kSize_);
      UTIL_CHECK(expKsq_.capacity() == kSize_);

      // Set up associated workspace fields slices
      fft().forwardTransform(qin, qk_);
      VecOp::mulEqV(qk_, expKsq_); // qk *= expKsq
      fft().inverseTransformUnsafe(qk_, qout); 
   }

   /*
   * Apply the field operator for the bead model.
   */
   template <int D>
   void Block<D>::stepFieldBead(RField<D>& q)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(expW_.capacity() == nx);
      UTIL_CHECK(q.capacity() == nx);

      VecOp::mulEqV(q, expW_); // q *= expW
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentrationThread(double prefactor)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isSolved());
      UTIL_CHECK(propagator(1).isSolved());
      UTIL_CHECK(cField().capacity() == nx);

      // Initialize cField to zero at all points
      VecOp::eqS(cField(), 0.0);

      Pscf::Rpg::Propagator<D> const & p0 = propagator(0);
      Pscf::Rpg::Propagator<D> const & p1 = propagator(1);

      addEqMulVVc(cField(), p0.q(0), p1.q(ns_ - 1), 1.0);
      addEqMulVVc(cField(), p0.q(ns_ - 1), p1.q(0), 1.0);

      for (int j = 1; j < ns_ - 1; j += 2) {
         // Odd indices
         addEqMulVVc(cField(), p0.q(j), p1.q(ns_ - 1 - j), 4.0);
      }
      for (int j = 2; j < ns_ - 2; j += 2) {
         // Even indices
         addEqMulVVc(cField(), p0.q(j), p1.q(ns_ - 1 - j), 2.0);
      }

      VecOp::mulEqS(cField(), (prefactor * ds_/3.0));
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentrationBead(double prefactor)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isSolved());
      UTIL_CHECK(propagator(1).isSolved());
      UTIL_CHECK(cField().capacity() == nx);

      // Initialize cField to zero at all points
      VecOp::eqS(cField(), 0.0);

      // References to forward and reverse propagators
      Pscf::Rpg::Propagator<D> const & p0 = propagator(0);
      Pscf::Rpg::Propagator<D> const & p1 = propagator(1);

      // Vertex 0 contribution (if owned by block)
      if (ownsVertex(0)) {
         addEqMulVVV(cField(), p0.q(0), p1.q(ns_ - 1), expWInv_);
      }

      // Internal beads
      for (int j = 1; j < ns_ - 1; ++j) {
         addEqMulVVV(cField(), p0.q(j), p1.q(ns_ - 1 - j), expWInv_);
      }

      // Vertex 0 contribution (if owned by block)
      if (ownsVertex(1)) {
         addEqMulVVV(cField(), p0.q(ns_-1), p1.q(0), expWInv_);
      }

      // Scale cField() by prefactor
      VecOp::mulEqS(cField(), prefactor);
   }


   /*
   * Average of a product of complementary propagator slices.
   *
   * This computes the spatial average of q0(r)*q1(r)*exp(+W(r))
   */
   template <int D>
   double
   Block<D>::averageProduct(RField<D> const& q0, RField<D> const& q1)
   {
      const int nx = mesh().size();
      UTIL_CHECK(q0.capacity() == nx);
      UTIL_CHECK(q1.capacity() == nx);

      double Q = Reduce::innerProduct(q0, q1);
      Q /= double(nx);
      return Q;
   }

   /*
   * Spatial integral of a product of complementary propagator slices.
   *
   * This computes the spatial average of q0(r)*q1(r)*exp(+W(r))
   */
   template <int D>
   double
   Block<D>::averageProductBead(RField<D> const& q0, RField<D> const& q1)
   {
      const int nx = mesh().size();
      UTIL_CHECK(q0.capacity() == nx);
      UTIL_CHECK(q1.capacity() == nx);
      if (qr_.isAllocated()) {
         UTIL_CHECK(qr_.capacity() == nx);
      } else {
         qr_.allocate(mesh().dimensions());
      }

      VecOp::mulVV(qr_, q0, q1);
      VecOp::mulEqV(qr_, expWInv_);
      double Q = Reduce::sum(qr_);
      Q /= double(nx);
      return Q;

   }

   /*
   * Compute stress contribution from this block.
   */
   template <int D>
   void Block<D>::computeStressThread(double prefactor)
   {
      int nx = mesh().size();

      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(kSize_ > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(nParams_ > 0);
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(propagator(0).isSolved());
      UTIL_CHECK(propagator(1).isSolved());

      // Calculate dKSq if necessary
      if (!waveListPtr_->hasdKSq()) {
         waveListPtr_->computedKSq();
      }

      // Workspace variables
      double dels, normal, increment;
      normal = 3.0*6.0;
      Pscf::Rpg::Propagator<D>& p0 = propagator(0);
      Pscf::Rpg::Propagator<D>& p1 = propagator(1);
      FSArray<double, 6> dQ;
      int i, j, n;
      RField<D> rTmp(kMeshDimensions_); // array of real values on kgrid

      // Initialize dQ and stress to 0
      stress_.clear();
      for (i = 0; i < nParams_; ++i) {
         dQ.append(0.0);
         stress_.append(0.0);
      }

      if (useBatchedFFT_) {

         UTIL_CHECK(fftBatchedAll_.isSetup());
         UTIL_CHECK(mesh().dimensions() == fftBatchedAll_.meshDimensions());
         UTIL_CHECK(!q0k_.isAllocated());
         UTIL_CHECK(!q1k_.isAllocated());
         // In this case, containers q0k_ and q1k_ will be associated
         // with slices of q0kBatched_, q1kBatched_

         // Computed batched FFT for propagators in both directions
         fftBatchedAll_.forwardTransform(p0.qAll(), q0kBatched_);
         fftBatchedAll_.forwardTransform(p1.qAll(), q1kBatched_);

      } else {

         if (!q0k_.isAllocated()) {
            q0k_.allocate(mesh().dimensions());
            q1k_.allocate(mesh().dimensions());
         }

      }

      // Main loop over contour points
      for (j = 0; j < ns_ ; ++j) {

         if (useBatchedFFT_) { 
            // Batched FFTs have already been computed
            // Associate q0k_, q1k_ with slices of q0kBatched_, q1kBatched_
            q0k_.associate(q0kBatched_, j * kSize_, mesh().dimensions());
            q1k_.associate(q1kBatched_, (ns_-1-j) * kSize_, 
                           mesh().dimensions());
         } else {
            // Compute Fourier transforms at contour grid point j
            UTIL_CHECK(fft().isSetup());
            fft().forwardTransform(p0.q(j), q0k_);
            fft().forwardTransform(p1.q(ns_-1-j), q1k_);
         }

         // Compute prefactor for Simpson's rule
         dels = ds_;
         if (j != 0 && j != ns_ - 1) {
            if (j % 2 == 0) {
               dels = dels*2.0;
            } else {
               dels = dels*4.0;
            }
         }

         // Increment stress contributions for all unit cell parameters
         for (n = 0; n < nParams_ ; ++n) {

            // Launch kernel to evaluate dQ at all wavevectors
            realMulVConjVV(rTmp, q0k_, q1k_, waveListPtr_->dKSq(n));

            // Get the sum of all elements
            increment = Reduce::sum(rTmp);
            increment *= kuhn() * kuhn() * dels / normal;
            dQ[n] -= increment;
         }

         if (useBatchedFFT_) {
            q0k_.dissociate();
            q1k_.dissociate();
         }

      } // end loop over contour points

      // Normalize total stress values
      for (i = 0; i < nParams_; ++i) {
         stress_[i] -= (dQ[i] * prefactor);
      }

   }

   /*
   * Compute stress contribution from this block, in bead model.
   */
   template <int D>
   void Block<D>::computeStressBead(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(PolymerModel::isBead());
      UTIL_CHECK(isAllocated_);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(kSize_ > 0);
      UTIL_CHECK(nParams_ > 0);
      UTIL_CHECK(mesh().dimensions() == fft().meshDimensions());
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isSolved());
      UTIL_CHECK(propagator(1).isSolved());

      // Calculate dKSq if necessary
      if (!waveListPtr_->hasdKSq()) {
         waveListPtr_->computedKSq();
      }

      int i, j, n;

      // Initialize dQ and stress to 0
      FSArray<double, 6> dQ;
      stress_.clear();
      for (i = 0; i < nParams_; ++i) {
         dQ.append(0.0);
         stress_.append(0.0);
      }

      // References to forward and reverse propagators
      Pscf::Rpg::Propagator<D>& p0 = propagator(0);
      Pscf::Rpg::Propagator<D>& p1 = propagator(1);

      if (useBatchedFFT_) {

         UTIL_CHECK(fftBatchedAll_.isSetup());
         UTIL_CHECK(mesh().dimensions() == fftBatchedAll_.meshDimensions());
         UTIL_CHECK(!q0k_.isAllocated());
         UTIL_CHECK(!q1k_.isAllocated());
         // In this case, containers q0k_ and q1k_ will be associated
         // with slices of q0kBatched_, q1kBatched_

         // Computed batched FFT for propagators in both directions
         fftBatchedAll_.forwardTransform(p0.qAll(), q0kBatched_);
         fftBatchedAll_.forwardTransform(p1.qAll(), q1kBatched_);

      } else {

         if (!q0k_.isAllocated()) {
            q0k_.allocate(mesh().dimensions());
            q1k_.allocate(mesh().dimensions());
         }

      }

      double increment = 0.0;
      double bSq = kuhn()*kuhn()/6.0;
      RField<D> rTmp(kMeshDimensions_); // array of real values on kgrid

      // Main loop over contour points
      for (j = 0; j < ns_ - 1; ++j) {

         if (useBatchedFFT_) { 
            // Batched FFTs have already been computed
            // Associate q0k_, q1k_ with slices of q0kBatched_, q1kBatched_
            q0k_.associate(q0kBatched_, j * kSize_, mesh().dimensions());
            q1k_.associate(q1kBatched_, (ns_- 2 -j) * kSize_, 
                           mesh().dimensions());
         } else {
            // Compute Fourier transforms at contour grid point j
            UTIL_CHECK(fft().isSetup());
            fft().forwardTransform(p0.q(j), q0k_);
            fft().forwardTransform(p1.q(ns_ - 2 - j), q1k_);
         }

         // Increment stress contributions for all unit cell parameters
         for (n = 0; n < nParams_ ; ++n) {

            // Launch kernel to evaluate dQ at all wavevectors
            realMulVConjVV(rTmp, q0k_, q1k_, waveListPtr_->dKSq(n));
            VecOp::mulEqV(rTmp, expKsq_);

            // Get the sum of all elements
            increment = Reduce::sum(rTmp);
            increment *= bSq;
            dQ[n] -= increment;
         }

         if (useBatchedFFT_) {
            q0k_.dissociate();
            q1k_.dissociate();
         }

      } // end loop over contour points

      // Normalize total stress values
      for (i = 0; i < nParams_; ++i) {
         stress_[i] -= dQ[i] * prefactor;
      }

   }

}
}
#endif
