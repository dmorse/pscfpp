#ifndef PSPG_BLOCK_TPP
#define PSPG_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/cuda/GpuResources.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/containers/FMatrix.h>       // member template
#include <util/containers/DArray.h>        // member template
#include <util/containers/FArray.h>        // member template
#include <sys/time.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Pscf::Prdc;

   // CUDA kernels (only used in this file)

   static __global__ 
   void mulDelKsq(cudaReal* result, const cudaComplex* q1,
                  const cudaComplex* q2, const cudaReal* delKsq,
                                    int paramN, int kSize, int rSize) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < kSize; i += nThreads) {
         #ifdef SINGLE_PRECISION
         result[i] =  cuCmulf( q1[i], 
                      cuConjf(q2[i])).x * delKsq[paramN * rSize + i];
         #else
         result[i] =  cuCmul( q1[i], 
                      cuConj(q2[i])).x * delKsq[paramN * rSize + i];
         #endif
      }
   }

   static __global__ 
   void pointwiseMulSameStart(const cudaReal* a, const cudaReal* expW,
                              const cudaReal* expW2,  
                              cudaReal* q1, cudaReal* q2, 
                              int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal input;
      for (int i = startID; i < size; i += nThreads) {
         input = a[i];
         q1[i] = expW[i] * input;
         q2[i] = expW2[i] * input;
      }
   }

   static __global__ 
   void pointwiseMulTwinned(const cudaReal* qr1, 
                            const cudaReal* qr2, 
                            const cudaReal* expW, 
                            cudaReal* q1, cudaReal* q2, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal scale;
      for (int i = startID; i < size; i += nThreads) {
         scale = expW[i];
         q1[i] = qr1[i] * scale;
         q2[i] = qr2[i] * scale;
      }
   }

   static __global__ 
   void scaleComplexTwinned(cudaComplex* qk1, cudaComplex* qk2, 
                            const cudaReal* expksq1, 
                            const cudaReal* expksq2, 
                            int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < size; i += nThreads) {
         qk1[i].x *= expksq1[i];
         qk1[i].y *= expksq1[i];
         qk2[i].x *= expksq2[i];
         qk2[i].y *= expksq2[i];
      }
   }

   static __global__ 
   void scaleComplex(cudaComplex* a, cudaReal* scale, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startID; i < size; i += nThreads) {
         a[i].x *= scale[i];
         a[i].y *= scale[i];
      }
   }

   static __global__ 
   void richardsonExpTwinned(cudaReal* qNew, const cudaReal* q1,
      const cudaReal* qr, const cudaReal* expW2, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal q2;
      for (int i = startID; i < size; i += nThreads) {
         q2 = qr[i] * expW2[i];
         qNew[i] = (4.0 * q2 - q1[i]) / 3.0;
      }
   }

   static __global__ 
   void multiplyScaleQQ(cudaReal* result,
                        const cudaReal* p1,
                        const cudaReal* p2,
                        double scale, int size) 
   {

      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;

      for(int i = startID; i < size; i += nThreads) {
         result[i] += scale * p1[i] * p2[i];
      }

   }

   // Block<D> member functions

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(0),
      kMeshDimensions_(0),
      ds_(0.0),
      ns_(0),
      temp_(0),
      isAllocated_(false),
      hasExpKsq_(false),
      expKsq_host(0),
      expKsq2_host(0)
   {
      propagator(0).setBlock(*this);
      propagator(1).setBlock(*this);
   }

   /*
   * Destructor.
   */
   template <int D>
   Block<D>::~Block()
   {      
      if (isAllocated_) {
         cudaFree(qkBatched_);
         cudaFree(qk2Batched_);
         cudaFree(d_temp_);
         delete[] temp_;
         delete[] expKsq_host;
         delete[] expKsq2_host;
         
      }
   }

   template <int D>
   void Block<D>::setDiscretization(double ds, Mesh<D> const & mesh)
   {
      UTIL_CHECK(mesh.size() > 1);
      UTIL_CHECK(ds > 0.0);
      UTIL_CHECK(!isAllocated_);

      // GPU Resources
      ThreadGrid::setThreadsLogical(mesh.size(),nBlocks_,nThreads_);

      // Set association to mesh
      meshPtr_ = &mesh;

      // Set contour length discretization for this block
      int tempNs;
      tempNs = floor( length()/(2.0 *ds) + 0.5 );
      if (tempNs == 0) {
         tempNs = 1;
      }
      ns_ = 2*tempNs + 1;

      ds_ = length()/double(ns_ - 1);

      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh.dimensions()[i];
         } else {
            kMeshDimensions_[i] = mesh.dimensions()[i]/2 + 1;
         }
      }

      kSize_ = 1;
      for(int i = 0; i < D; ++i) {
         kSize_ *= kMeshDimensions_[i];
      }

      // Allocate work arrays
      expKsq_.allocate(kMeshDimensions_);
      expKsq2_.allocate(kMeshDimensions_);
      expW_.allocate(mesh.dimensions());
      expW2_.allocate(mesh.dimensions());
      qr_.allocate(mesh.dimensions());
      qr2_.allocate(mesh.dimensions());
      qk_.allocate(mesh.dimensions());
      qk2_.allocate(mesh.dimensions());
      q1_.allocate(mesh.dimensions());
      q2_.allocate(mesh.dimensions());

      propagator(0).allocate(ns_, mesh);
      propagator(1).allocate(ns_, mesh);
      gpuErrchk( cudaMalloc((void**)&qkBatched_, 
                             ns_ * kSize_ * sizeof(cudaComplex))
               );
      gpuErrchk( cudaMalloc((void**)&qk2Batched_, 
                             ns_ * kSize_ * sizeof(cudaComplex))
               );
      cField().allocate(mesh.dimensions());

      gpuErrchk(cudaMalloc((void**)&d_temp_, nBlocks_ * sizeof(cudaReal)));
      temp_ = new cudaReal[nBlocks_];

      expKsq_host = new cudaReal[kSize_];
      expKsq2_host = new cudaReal[kSize_];

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
   void Block<D>::setupUnitCell(const UnitCell<D>& unitCell, 
                                const WaveList<D>& wavelist)
   {
      nParams_ = unitCell.nParameter();

      // store pointer to unit cell and wavelist
      unitCellPtr_ = &unitCell;
      waveListPtr_ = &wavelist;

      hasExpKsq_ = false;
   }
   
   /*
   * Setup data that depend on the unit cell parameters.
   */
   template <int D>
   void Block<D>::setupUnitCell(const UnitCell<D>& unitCell)
   {
      nParams_ = unitCell.nParameter();

      // store pointer to unit cell and wavelist
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

      // Setup expKsq values on Host then transfer to device
      int kSize = 1;
      for(int i = 0; i < D; ++i) {
         kSize *= kMeshDimensions_[i];
      }

      int i;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank();
         if (waveListPtr_ == 0){
            G = iter.position();
            Gmin = shiftToMinimum(G, mesh().dimensions(), unitCell());
            Gsq = unitCell().ksq(Gmin);
         } else{
            Gsq = unitCell().ksq(wavelist().minImage(iter.rank()));
         }
         expKsq_host[i] = exp(Gsq*factor);
         expKsq2_host[i] = exp(Gsq*factor / 2);
      }

      cudaMemcpy(expKsq_.cField(), expKsq_host, 
                 kSize * sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(expKsq2_.cField(), expKsq2_host, 
                 kSize * sizeof(cudaReal), cudaMemcpyHostToDevice);

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
      assignExp<<<nBlocks_, nThreads_>>>(expW_.cField(), w.cField(), 
                                         (double)0.5* ds_, nx);
      assignExp<<<nBlocks_, nThreads_>>>(expW2_.cField(), w.cField(), 
                                         (double)0.25 * ds_, nx);

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
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == nx)

      // Initialize cField to zero at all points
      assignUniformReal<<<nBlocks_, nThreads_>>>
              (cField().cField(), 0.0, nx);

      Pscf::Pspg::Propagator<D> const & p0 = propagator(0);
      Pscf::Pspg::Propagator<D> const & p1 = propagator(1);

      multiplyScaleQQ<<<nBlocks_, nThreads_>>>
              (cField().cField(), p0.q(0), p1.q(ns_ - 1), 1.0, nx);
      multiplyScaleQQ<<<nBlocks_, nThreads_>>>
              (cField().cField(), p0.q(ns_-1), p1.q(0), 1.0, nx);
      for (int j = 1; j < ns_ - 1; j += 2) {
         // Odd indices
         multiplyScaleQQ<<<nBlocks_, nThreads_>>>
                (cField().cField(), p0.q(j), p1.q(ns_ - 1 - j), 4.0, nx);
      }
      for (int j = 2; j < ns_ - 2; j += 2) {
          // Even indices
          multiplyScaleQQ<<<nBlocks_, nThreads_>>>
                (cField().cField(), p0.q(j), p1.q(ns_ - 1 - j), 2.0, nx);
      }

      scaleReal<<<nBlocks_, nThreads_>>>
          (cField().cField(), (prefactor * ds_/3.0), nx);

   }

   template <int D>
   void Block<D>::setupFFT() {
      if (!fft_.isSetup()) {
         fft_.setup(qr_, qk_);
         fftBatched_.setup(mesh().dimensions(), kMeshDimensions_, ns_);
      }
   }

   /*
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::step(const cudaReal* q, cudaReal* qNew)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasExpKsq_);

      // Check real-space mesh sizes
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);

      // Fourier-space mesh sizes
      int nk = qk_.capacity();
      UTIL_CHECK(expKsq_.capacity() == nk);

      // Apply pseudo-spectral algorithm

      pointwiseMulSameStart<<<nBlocks_, nThreads_>>>
                           (q, expW_.cField(), expW2_.cField(), 
                            qr_.cField(), qr2_.cField(), nx);
      fft_.forwardTransform(qr_, qk_);
      fft_.forwardTransform(qr2_, qk2_);
      scaleComplexTwinned<<<nBlocks_, nThreads_>>>
                         (qk_.cField(), qk2_.cField(), 
                          expKsq_.cField(), expKsq2_.cField(), nk);
      fft_.inverseTransform(qk_, qr_);
      fft_.inverseTransform(qk2_, q2_);
      pointwiseMulTwinned<<<nBlocks_, nThreads_>>>
                         (qr_.cField(), q2_.cField(), expW_.cField(), 
                          q1_.cField(), qr_.cField(), nx);
      fft_.forwardTransform(qr_, qk_);
      scaleComplex<<<nBlocks_, nThreads_>>>(qk_.cField(), expKsq2_.cField(), nk);
      fft_.inverseTransform(qk_, qr_);
      richardsonExpTwinned<<<nBlocks_, nThreads_>>>(qNew, q1_.cField(),
                           qr_.cField(), expW2_.cField(), nx);
      //remove the use of q2


   }


   /*
   * Compute stress contribution from this block. 
   */
   template <int D>
   void Block<D>::computeStress(WaveList<D> const & wavelist, 
                                double prefactor)
   {
      int nx = mesh().size();

      // Preconditions
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      double dels, normal, increment;
      normal = 3.0*6.0;

      FArray<double, 6> dQ;
      int i;
      for (i = 0; i < 6; ++i) {
         dQ [i] = 0.0;
         stress_[i] = 0.0;
      }

      Pscf::Pspg::Propagator<D> const & p0 = propagator(0);
      Pscf::Pspg::Propagator<D> const & p1 = propagator(1);

      fftBatched_.forwardTransform(p0.head(), qkBatched_, ns_);
      fftBatched_.forwardTransform(p1.head(), qk2Batched_, ns_);
      cudaMemset(qr2_.cField(), 0, mesh().size() * sizeof(cudaReal));

      for (int j = 0; j < ns_ ; ++j) {

         dels = ds_;
         if (j != 0 && j != ns_ - 1) {
            if (j % 2 == 0) {
               dels = dels*2.0;
            } else {
               dels = dels*4.0;
            }
         }

         for (int n = 0; n < nParams_ ; ++n) {
            mulDelKsq<<<nBlocks_, nThreads_ >>>
                (qr2_.cField(), 
                 qkBatched_ + (j*kSize_), 
                 qk2Batched_ + (kSize_ * (ns_ -1 -j)),
                 wavelist.dkSq(), n , kSize_, nx);

            increment = gpuSum(qr2_.cField(), mesh().size());
            increment = (increment * kuhn() * kuhn() * dels)/normal;
            dQ [n] = dQ[n]-increment;
         }

      }

      // Normalize
      for (i = 0; i < nParams_; ++i) {
         stress_[i] = stress_[i] - (dQ[i] * prefactor);
      }
   }

}
}
#endif
