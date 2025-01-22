#ifndef RPG_BLOCK_TPP
#define RPG_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/cuda/GpuResources.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

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
      * Performs qNew = (4 * (qr2 * expW2) - qr) / 3 elementwise, CUDA kernel
      */
      __global__ void _richardsonEx(cudaReal* qNew, cudaReal const * qr,
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
      __global__ void _addEqMulVVc(cudaReal* a, cudaReal const * b,
                                   cudaReal const * c, cudaReal const d, 
                                   const int n) 
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for(int i = startID; i < n; i += nThreads) {
            a[i] += b[i] * c[i] * d;
         }
      }

   }

   // CUDA kernel wrappers:

   /*
   * Element-wise calculation of a = real(b * conj(c) * d), kernel wrapper
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
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _realMulVConjVV<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), 
                                             c.cArray(), d.cArray(), n);
   }

   /*
   * Performs qNew = (4 * (qr2 * expW2) - qr) / 3 elementwise, kernel wrapper
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
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _richardsonEx<<<nBlocks, nThreads>>>(qNew.cArray(), qr.cArray(), 
                                           qr2.cArray(), expW2.cArray(), n);
   }

   /*
   * Performs a[i] += b[i] * c[i] * d, kernel wrapper
   */
   void addEqMulVVc(DeviceArray<cudaReal>& a, DeviceArray<cudaReal> const & b,
                    DeviceArray<cudaReal> const & c, cudaReal const d) 
   {
      int n = a.capacity();
      UTIL_CHECK(b.capacity() >= n);
      UTIL_CHECK(c.capacity() >= n);
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      // Launch kernel
      _addEqMulVVc<<<nBlocks, nThreads>>>(a.cArray(), b.cArray(), 
                                          c.cArray(), d, n);
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

      // store pointers to associated objects
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      unitCellPtr_ = &cell;
      waveListPtr_ = &wavelist;

      // Compute Fourier space kMeshDimensions_ and kSize_
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh.dimensions()[i];
         } else {
            kMeshDimensions_[i] = mesh.dimensions()[i]/2 + 1;
         }
         kSize_ *= kMeshDimensions_[i];
      }

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

      // Set contour length discretization for this block
      dsTarget_ = ds;
      int tempNs;
      tempNs = floor(length() / (2.0 * ds) + 0.5);
      if (tempNs == 0) {
         tempNs = 1; // ensure at least 3 contour steps per chain
      }
      ns_ = 2*tempNs + 1;

      ds_ = length()/double(ns_ - 1);

      // Setup fftBatched objects
      UTIL_CHECK(!fftBatchedPair_.isSetup());
      fftBatchedPair_.setup(mesh().dimensions(), 2);
      if (useBatchedFFT_) {
         UTIL_CHECK(!fftBatchedAll_.isSetup());
         fftBatchedAll_.setup(mesh().dimensions(), ns_);
      }

      // Allocate work arrays
      expKsq_.allocate(kMeshDimensions_);
      expKsq2_.allocate(kMeshDimensions_);
      expW_.allocate(mesh().dimensions());
      expW2_.allocate(mesh().dimensions());
      qrPair_.allocate(2 * mesh().size());
      qkPair_.allocate(2 * kSize_);
      q1_.allocate(mesh().dimensions());
      q2_.allocate(mesh().dimensions());

      propagator(0).allocate(ns_, mesh());
      propagator(1).allocate(ns_, mesh());
      
      cField().allocate(mesh().dimensions());

      if (useBatchedFFT_) {
         qkBatched_.allocate(ns_ * kSize_);
         qk2Batched_.allocate(ns_ * kSize_);
      }

      expKsq_h_.allocate(kSize_);
      expKsq2_h_.allocate(kSize_);

      isAllocated_ = true;
      hasExpKsq_ = false;
   }

   /*
   * Upon changing lattice parameters, update this object.
   */
   template <int D>
   void Block<D>::updateUnitCell()
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
      BlockDescriptor::setLength(newLength);
      
      if (isAllocated_) { // if allocate() has already been called
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
               qkBatched_.deallocate();
               qk2Batched_.deallocate();
               qkBatched_.allocate(ns_ * kSize_);
               qk2Batched_.allocate(ns_ * kSize_);
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
      
      double factor = -1.0*kuhn()*kuhn()*ds_/6.0;

      // Calculate expKsq values on device
      VecOp::expVc(expKsq_, waveListPtr_->kSq(), factor);
      VecOp::expVc(expKsq2_, waveListPtr_->kSq(), factor / 2.0);

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
      VecOp::expVc(expW_, w, -0.5 * ds_);
      VecOp::expVc(expW2_, w, -0.25 * ds_);

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
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::step(RField<D> const & q, RField<D>& qNew)
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

      // Set up associated workspace fields
      RField<D> qr, qr2;
      RFieldDft<D> qk, qk2;
      qr.associate(qrPair_, 0, mesh().dimensions());
      qr2.associate(qrPair_, nx, mesh().dimensions());
      qk.associate(qkPair_, 0, mesh().dimensions());
      qk2.associate(qkPair_, kSize_, mesh().dimensions());

      // Apply pseudo-spectral algorithm
      VecOp::mulVVPair(qr, qr2, expW_, expW2_, q); // qr = expW*q, qr2 = expW2*q
      fftBatchedPair_.forwardTransform(qrPair_, qkPair_); // to Fourier space
      VecOp::mulEqV(qk, expKsq_); // qk *= expKsq
      VecOp::mulEqV(qk2, expKsq2_); // qk2 *= expKsq2
      fftBatchedPair_.inverseTransformUnsafe(qkPair_, qrPair_); // to real space
      VecOp::mulEqVPair(qr, qr2, expW_); // qr *= expW, qr2 *= expW
      fft().forwardTransform(qr2, qk2); // to Fourier space, only qr2
      VecOp::mulEqV(qk2, expKsq2_); // qk2 *= expKsq2
      fft().inverseTransformUnsafe(qk2, qr2); // to real space, only qr2
      richardsonEx(qNew, qr, qr2, expW2_); // qNew=(4*(qr2*expW2)-qr)/3
   }

   /*
   * Compute stress contribution from this block. 
   */
   template <int D>
   void Block<D>::computeStress(double prefactor)
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
      RFieldDft<D> qk, qk2;

      // Initialize dQ and stress to 0
      stress_.clear();
      for (i = 0; i < nParams_; ++i) {
         dQ.append(0.0);
         stress_.append(0.0);
      }

      if (useBatchedFFT_) {
         // Get q at all contour points in Fourier space for both propagators
         UTIL_CHECK(fftBatchedAll_.isSetup());
         UTIL_CHECK(mesh().dimensions() == fftBatchedAll_.meshDimensions());
         fftBatchedAll_.forwardTransform(p0.qAll(), qkBatched_);
         fftBatchedAll_.forwardTransform(p1.qAll(), qk2Batched_);
      } else {
         // Allocate qk and qk2 to store results of individual FFTs
         qk.allocate(mesh().dimensions());
         qk2.allocate(mesh().dimensions());
      }

      // Main loop over contour points
      for (j = 0; j < ns_ ; ++j) {

         if (useBatchedFFT_) { // FFTs have already been calculated
            // Associate qk and qk2 with sections of qkBatched_ and qk2Batched_
            qk.associate(qkBatched_, j * kSize_, mesh().dimensions());
            qk2.associate(qk2Batched_, (ns_-1-j) * kSize_, mesh().dimensions());
         } else {
            // Get q at contour point j in Fourier space for both propagators
            UTIL_CHECK(fft().isSetup());
            fft().forwardTransform(p0.q(j), qk);
            fft().forwardTransform(p1.q(ns_-1-j), qk2);
         }

         dels = ds_;
         if (j != 0 && j != ns_ - 1) {
            if (j % 2 == 0) {
               dels = dels*2.0;
            } else {
               dels = dels*4.0;
            }
         }

         for (n = 0; n < nParams_ ; ++n) {
            // Launch kernel to evaluate dQ for each basis function
            realMulVConjVV(rTmp, qk, qk2, waveListPtr_->dKSq(n));

            // Get the sum of all elements
            increment = Reduce::sum(rTmp);
            increment *= kuhn() * kuhn() * dels / normal;
            dQ[n] -= increment;
         }

         if (useBatchedFFT_) {
            qk.dissociate();
            qk2.dissociate();
         }
      }

      // Normalize
      for (i = 0; i < nParams_; ++i) {
         stress_[i] -= (dQ[i] * prefactor);
      }
   }

}
}
#endif
