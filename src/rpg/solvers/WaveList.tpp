#ifndef RPG_WAVE_LIST_TPP
#define RPG_WAVE_LIST_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"
#include <pscf/cuda/GpuResources.h>
#include <pscf/mesh/MeshIterator.h>

namespace Pscf {
namespace Rpg
{

   // CUDA kernels: 
   // (defined in anonymous namespace, used only in this file)

   namespace {

      /*
      * Compute minimum images of each wavevector, template declaration.
      */
      template <int D>
      __global__ void _computeMinimumImages(int* minImages, cudaReal* kSq,
                                            cudaReal const * kBasis, 
                                            int const * meshDims, 
                                            int const kSize);
      
      /*
      * Compute minimum images of each wavevector, 1D explicit specialization.
      *
      * Kernel must be launched with kSize threads. (One thread calculates 
      * one minimum image.)
      * 
      * When launched, this kernel requires no dynamic memory allocation.
      */
      template <>
      __global__ void _computeMinimumImages<1>(int* minImages, cudaReal* kSq,
                                               cudaReal const * kBasis, 
                                               int const * meshDims, 
                                               int const kSize)
      {
         unsigned int nThreads = blockDim.x * gridDim.x;
         unsigned int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load kBasis and meshDims
         cudaReal kBasis_ = kBasis[0];
         int meshDims_ = meshDims[0];

         for (int i = startID; i < kSize; i += nThreads) {
            int img = minImages[i];

            while (img > (meshDims_>>1)) { // note: x>>1 is same as x/2
               img -= meshDims_;
            }
            while (img < -1*(meshDims_>>1)) { 
               img += meshDims_;
            }

            minImages[i] = img;
            kSq[i] = img * img * kBasis_ * kBasis_;
         }
      }

      /*
      * Compute minimum images of each wavevector, 2D explicit specialization.
      *
      * This kernel should be launched with >=32*kSize threads.
      *
      * 32 threads are used to calculate one minimum image. In the CPU code,
      * we compare 25 different images of a wave to determine the minimum 
      * image, and here we choose 32 instead because it is a power of 2, 
      * allowing us to use a reduction algorithm to compare the images.
      * 
      * When launched, this kernel requires dynamic memory allocation of
      * (2*sizeof(int) + sizeof(cudaReal)) * nThreadsPerBlock.
      */
      template <>
      __global__ void _computeMinimumImages<2>(int* minImages, cudaReal* kSq,
                                               cudaReal const * kBasis, 
                                               int const * meshDims, 
                                               int const kSize)
      {
         unsigned int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Determine which lattice parameter and which image to evaluate
         unsigned int paramID = (startID >> 5); // equivalent to startID / 32
         unsigned int imageID = startID - (paramID * 32); // value in [0, 31]
         
         // Determine the image that will be evaluated by this thread
         // (uses integer division, not ideal for speed)
         int s0, s1;
         s0 = imageID / 5;
         s1 = imageID - (s0 * 5);
         s0 -= 2; // shift values to go from -2 to 4, not 0 to 6
         s1 -= 2; // shift values to go from -2 to 2, not 0 to 4

         // Set epsilon
         #ifdef SINGLE_PRECISION
         cudaReal epsilon = 1.0E-4;
         #else
         cudaReal epsilon = 1.0E-8;
         #endif

         // Declare dynamically allocated arrays in shared memory
         // (allocated space is shared between cudaReal array and int array)
         extern __shared__ cudaReal kSqVals_[];
         int* images_ = (int*)&kSqVals_[blockDim.x];

         if (paramID < kSize) { // only evaluate if on the k-grid

            // Load image from global memory and shift based on s0 and s1
            images_[2*threadIdx.x] = minImages[paramID] + (s0 * meshDims[0]);
            images_[2*threadIdx.x+1] = minImages[kSize + paramID] + 
                                       (s1 * meshDims[1]);

            // Calculate kSq for this wave
            cudaReal kVec0 = (kBasis[0] * images_[2*threadIdx.x]) + 
                             (kBasis[2] * images_[2*threadIdx.x+1]);
            cudaReal kVec1 = (kBasis[1] * images_[2*threadIdx.x]) + 
                             (kBasis[3] * images_[2*threadIdx.x+1]);

            kSqVals_[threadIdx.x] = (kVec0*kVec0) + (kVec1*kVec1);
         }
         __syncthreads(); // wait for all threads to finish

         // Perform a parallel reduction on 32 threads to find min image
         // (note: stride >>= 1 is equivalent to stride /= 2)
         for (int stride = 16; stride > 0; stride >>= 1) {
            if (paramID < kSize) { // only evaluate if on the k-grid
               if (imageID < stride) {
                  bool swap = false;
                  if (kSqVals_[threadIdx.x+stride] < 
                                       (kSqVals_[threadIdx.x] - epsilon)) {
                     swap = true;
                  } else if (kSqVals_[threadIdx.x+stride] < 
                                       (kSqVals_[threadIdx.x] + epsilon)) {
                     // kSq values effectively equal. 
                     // Determine whether to swap based on hkl indices
                     for (int i = 0; i < 2; ++i) {
                        if (images_[2*threadIdx.x+i] > 
                                       images_[2*(threadIdx.x+stride)+i]) {
                           break;
                        } else
                        if (images_[2*threadIdx.x+i] < 
                                       images_[2*(threadIdx.x+stride)+i]) {
                           swap = true;
                           break;
                        }
                     }
                  }
                  if (swap) {
                     images_[2*threadIdx.x] = 
                                       images_[2*(threadIdx.x+stride)];
                     images_[2*threadIdx.x+1] = 
                                       images_[2*(threadIdx.x+stride)+1];
                     kSqVals_[threadIdx.x] = kSqVals_[threadIdx.x+stride];
                  }
               }
            }
            // note: no __syncthreads() needed here because this reduction 
            // occurswithin 1 warp, so all threads are already synced)
         }

         // At this point, for any thread with imageID == 0, the corresponding
         // entries in kSqVals_ and images_ should contain the kSq value and 
         // the minimum image, respectively. Store values and exit
         if ((imageID == 0) && (paramID < kSize)) {
            kSq[paramID] = kSqVals_[threadIdx.x];
            minImages[paramID] = images_[2*threadIdx.x];
            minImages[paramID + kSize] = images_[2*threadIdx.x+1];
         }
      }

      /*
      * Compute minimum images of each wavevector, 3D explicit specialization.
      *
      * This kernel should be launched with >=128*kSize threads.
      * 
      * 128 threads are used to calculate one minimum image. In the CPU code,
      * we compare 125 different images of a wave to determine the minimum 
      * image, and here we choose 128 instead because it is a power of 2, 
      * allowing us to use a reduction algorithm to compare the images.
      * 
      * When launched, this kernel requires dynamic memory allocation of
      * (3*sizeof(int) + sizeof(cudaReal)) * nThreadsPerBlock.
      */
      template <>
      __global__ void _computeMinimumImages<3>(int* minImages, cudaReal* kSq,
                                               cudaReal const * kBasis, 
                                               int const * meshDims, 
                                               int const kSize)
      {
         unsigned int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Determine which lattice parameter and which image to evaluate
         unsigned int paramID = (startID >> 7); // equivalent to startID / 128
         unsigned int imageID = startID - (paramID * 128);

         // Determine the image that will be evaluated by this thread
         // (uses integer division, not ideal for speed)
         int s0, s1, s2;
         s0 = imageID / 25;
         s1 = (imageID - (s0 * 25)) / 5;
         s2 = imageID - (s0 * 25) - (s1 * 5);
         s0 -= 2; // shift values to go from -2 to 2
         s1 -= 2;
         s2 -= 2;

         // Set epsilon
         #ifdef SINGLE_PRECISION
         cudaReal epsilon = 1.0E-4;
         #else
         cudaReal epsilon = 1.0E-8;
         #endif

         // Declare dynamically allocated arrays in shared memory
         // (allocated space is shared between cudaReal array and int array)
         extern __shared__ cudaReal kSqVals_[];
         int* images_ = (int*)&kSqVals_[blockDim.x];

         if (paramID < kSize) { // only evaluate if on the k-grid

            // Load image data from global memory
            images_[3*threadIdx.x] = minImages[paramID] + (s0 * meshDims[0]);
            images_[3*threadIdx.x+1] = minImages[kSize + paramID] + 
                                       (s1 * meshDims[1]);
            images_[3*threadIdx.x+2] = minImages[kSize + kSize + paramID] + 
                                       (s2 * meshDims[2]);

            // Calculate kSq for this wave
            cudaReal kVec0(0.0), kVec1(0.0), kVec2(0.0);
            for (int k = 0; k < 3; ++k) {
               kVec0 += kBasis[3*k]   * images_[3*threadIdx.x + k];
               kVec1 += kBasis[3*k+1] * images_[3*threadIdx.x + k];
               kVec2 += kBasis[3*k+2] * images_[3*threadIdx.x + k];
            }

            kSqVals_[threadIdx.x] = (kVec0*kVec0) + (kVec1*kVec1) + 
                                    (kVec2*kVec2);
         }
         __syncthreads(); // wait for all threads to finish

         // Perform a parallel reduction on 128 threads to find min image
         // (note: stride >>= 1 is equivalent to stride /= 2)
         for (int stride = 64; stride > 0; stride >>= 1) {
            if (paramID < kSize) { // only evaluate if on the k-grid
               if (imageID < stride) {
                  bool swap = false;
                  if (kSqVals_[threadIdx.x+stride] < 
                                          (kSqVals_[threadIdx.x] - epsilon)) {
                     swap = true;
                  } else if (kSqVals_[threadIdx.x+stride] < 
                                          (kSqVals_[threadIdx.x] + epsilon)) {
                     // kSq values effectively equal. 
                     // Determine whether to swap based on hkl indices
                     for (int i = 0; i < 3; ++i) {
                        if (images_[3*threadIdx.x+i] > 
                                       images_[3*(threadIdx.x+stride)+i]) {
                           break;
                        } else
                        if (images_[3*threadIdx.x+i] < 
                                       images_[3*(threadIdx.x+stride)+i]) {
                           swap = true;
                           break;
                        }
                     }
                  }
                  if (swap) {
                     images_[3*threadIdx.x] = 
                                          images_[3*(threadIdx.x+stride)];
                     images_[3*threadIdx.x+1] = 
                                          images_[3*(threadIdx.x+stride)+1];
                     images_[3*threadIdx.x+2] = 
                                          images_[3*(threadIdx.x+stride)+2];
                     kSqVals_[threadIdx.x] = kSqVals_[threadIdx.x+stride];
                  }
               }
            }
            __syncthreads(); // wait for all threads to finish
         }

         // At this point, for any thread with imageID == 0, the corresponding
         // entries in kSqVals_ and images_ should contain the kSq value and 
         // the minimum image, respectively. Store values and exit
         if ((imageID == 0) && (paramID < kSize)) {
            kSq[paramID] = kSqVals_[threadIdx.x];
            minImages[paramID] = images_[3*threadIdx.x];
            minImages[paramID + kSize] = images_[3*threadIdx.x+1];
            minImages[paramID + kSize + kSize] = images_[3*threadIdx.x+2];
         }
      }

      /*
      * Compute the kSq array on the GPU.
      */
      template <int D>
      __global__ void _computeKSq(cudaReal* kSq, int const * waveBz,
                                  cudaReal const * kBasis,
                                  int const nParams, int const kSize) 
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load kBasis into shared memory for fast access
         __shared__ cudaReal kBasis_[D*D];
         if (threadIdx.x < D * D) {
            kBasis_[threadIdx.x] = kBasis[threadIdx.x];
         }
         __syncthreads(); // wait for all threads to finish

         // Variables to be used in the loop
         int i, j, k, waveBz_;
         cudaReal kVec[D], kSqVal;
         // Note: usually local arrays are very slow in a CUDA kernel, but 
         // the compiler should ideally be able to unroll the loops below 
         // and store the kVec array in the register, because the full
         // structure of the loops below is known at compile time.
         
         // Loop through array
         for (i = startID; i < kSize; i += nThreads) {

            // initialize kVec to 0
            for (j = 0; j < D; ++j) {
               kVec[j] = 0.0;
            }

            // Calculate kVec
            for (j = 0; j < D; ++j) {
               waveBz_ = waveBz[i + (j * kSize)];
               for (k = 0; k < D; ++k) {
                  kVec[k] += kBasis_[k + (D*j)] * waveBz_;
               }
            }

            // Compute kSq
            kSqVal = 0.0;
            for (j = 0; j < D; ++j) {
               kSqVal += kVec[j] * kVec[j];
            }

            // Store value in global memory
            kSq[i] = kSqVal;

         } // kSize
      }

      /*
      * Compute the dKSq array on the GPU
      */
      template <int D>
      __global__ void _computedKSq(cudaReal* dKSq, int const * waveBz,
                                   cudaReal const * dkkBasis,
                                   bool const * implicitInverse, 
                                   int const nParams, int const kSize) 
      {
         // Size of dKSq is kSize * nParams
         // Each thread does nParams calculations
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load dkkBasis into shared memory for fast access
         // (max size of dkkBasis is 54 elements for triclinic unit cell)
         extern __shared__ cudaReal dkkBasis_[];
         if (threadIdx.x < nParams * D * D) {
            dkkBasis_[threadIdx.x] = dkkBasis[threadIdx.x];
         }
         __syncthreads(); // wait for all threads to finish

         // Variables to be used in the loop
         int param, i, j, k, waveBz_[D];
         cudaReal dKSqVal;
         
         // Loop through array
         for (param = 0; param < nParams; ++param) {
            for (i = startID; i < kSize; i += nThreads) {

               // initialize to 0
               dKSqVal = 0.0;

               // Load waveBz to local memory
               for (j = 0; j < D; ++j) {
                  waveBz_[j] = waveBz[i + (j * kSize)];
               }

               // Compute dKSq
               for (j = 0; j < D; ++j) {
                  for (k = 0; k < D; ++k) {
                     dKSqVal += waveBz_[j] * waveBz_[k]
                                * dkkBasis_[k + (j * D) + (param * D * D)];
                  } // D
               } // D

               if (implicitInverse[i]) { // if element i's inverse is implicit
                  dKSqVal *= 2;
               }

               dKSq[(param * kSize) + i] = dKSqVal;

            } // kSize
         } // nParams
      }
   }

   template <int D>
   WaveList<D>::WaveList()
    : kSize_(0),
      isAllocated_(false),
      hasMinimumImages_(false),
      hasKSq_(false),
      hasdKSq_(false),
      unitCellPtr_(nullptr),
      meshPtr_(nullptr)
   {}

   template <int D>
   WaveList<D>::~WaveList() 
   {}

   template <int D>
   void WaveList<D>::allocate(Mesh<D> const & m, UnitCell<D> const & c) 
   {
      UTIL_CHECK(m.size() > 0);
      UTIL_CHECK(c.nParameter() > 0);
      UTIL_CHECK(!isAllocated_);

      // Create permanent associations with mesh and unit cell
      unitCellPtr_ = &c;
      meshPtr_ = &m;

      int nParams = unitCell().nParameter();

      // Compute DFT mesh size kSize_ and dimensions kMeshDimensions_
      kSize_ = 1;
      for(int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh().dimension(i);
         } else {
            kMeshDimensions_[i] = mesh().dimension(i) / 2 + 1;
         }
         kSize_ *= kMeshDimensions_[i];
      }

      minImages_.allocate(kSize_ * D);
      kSq_.allocate(kMeshDimensions_);
      dKSq_.allocate(kSize_ * nParams);
      dKSqSlices_.allocate(nParams);
      for (int i = 0; i < nParams; i++) {
         dKSqSlices_[i].associate(dKSq_, i*kSize_, kMeshDimensions_);
      }

      // Set up implicitInverse_ array (only depends on mesh dimensions)
      implicitInverse_.allocate(kSize_);
      MeshIterator<D> kItr(kMeshDimensions_);
      HostDArray<bool> implicitInverse_h(kSize_);
      int inverseId;
      for (kItr.begin(); !kItr.atEnd(); ++kItr) {
         if (kItr.position(D-1) == 0) {
            inverseId = 0;
         } else {
            inverseId = mesh().dimension(D-1) - kItr.position(D-1);
         }
         if (inverseId > kMeshDimensions_[D-1]) {
            implicitInverse_h[kItr.rank()] = true;
         } else {
            implicitInverse_h[kItr.rank()] = false;
         }
      }
      implicitInverse_ = implicitInverse_h; // transfer to device memory

      isAllocated_ = true;
   }

   template <int D>
   void WaveList<D>::updateUnitCell()
   {
      if (hasVariableAngle()) {
         hasMinimumImages_ = false;
      }
      hasKSq_ = false;
      hasdKSq_ = false;
   }

   template <int D>
   void WaveList<D>::computeMinimumImages() 
   {
      if (hasMinimumImages_) return; // min images already calculated

      // Precondition
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(minImages_.capacity() == kSize_ * D);

      // Set initial array of images to contain the k-grid points
      HostDArray<int> imagesTmp(D*kSize_);
      MeshIterator<D> kItr(kMeshDimensions_);
      for (int i = 0; i < D; i++) {
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            imagesTmp[kItr.rank() + (i*kSize_)] = kItr.position(i);
         }
      }
      minImages_ = imagesTmp; // copy to device

      // Get kBasis and meshDims and store on device
      HostDArray<cudaReal> kBasis_h(D*D);
      HostDArray<int> meshDims_h(D);
      DeviceArray<cudaReal> kBasis(D*D);
      DeviceArray<int> meshDims(D);
      int idx = 0;
      for(int j = 0; j < D; ++j) {
         for(int k = 0; k < D; ++k) {
            kBasis_h[idx] = unitCell().kBasis(j)[k];
            idx++;
         }
         meshDims_h[j] = mesh().dimension(j);
      }
      kBasis = kBasis_h;
      meshDims = meshDims_h;

      // Set number of threads per gridpoint (depends on D)
      int threadsPerGP;
      if (D == 3) {
         threadsPerGP = 128;
      } else if (D == 2) {
         threadsPerGP = 32;
      } else if (D == 1) {
         threadsPerGP = 1;
      }

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(kSize_*threadsPerGP, nBlocks, nThreads);

      if ((D == 3) && (nThreads < 128)) {
         // Thread blocks too small. Manually set nThreads to 128
         ThreadGrid::setThreadsPerBlock(128);
         ThreadGrid::setThreadsLogical(kSize_*threadsPerGP, nBlocks, nThreads);

         // If the above was successful, print warning
         Log::file() << "Warning: " 
                     << "nThreads too small for computeMinimumImages.\n"
                     << "Setting nThreads equal to 128." << std::endl;
      }

      // Launch kernel
      size_t sz = (D * sizeof(int) + sizeof(cudaReal)) * nThreads;
      _computeMinimumImages<D><<<nBlocks, nThreads, sz>>>
         (minImages_.cArray(), kSq_.cArray(), kBasis.cArray(), 
          meshDims.cArray(), kSize_);
      
      hasMinimumImages_ = true;
      hasKSq_ = true;
   }

   template <int D>
   void WaveList<D>::computeKSq() 
   {
      if (hasKSq_) return; // kSq already calculated

      if (!hasMinimumImages_) {
         computeMinimumImages(); // compute both min images and kSq
         return;
      }

      // If this point is reached, calculate kSq using _computeKSq kernel

      // Precondition
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());

      // Get kBasis and store on device
      HostDArray<cudaReal> kBasis_h(D*D);
      DeviceArray<cudaReal> kBasis(D*D);
      int idx = 0;
      for(int j = 0; j < D; ++j) {
         for(int k = 0; k < D; ++k) {
            kBasis_h[idx] = unitCell().kBasis(j)[k];
            idx++;
         }
      }
      kBasis = kBasis_h;

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(kSize_, nBlocks, nThreads);
      
      // Launch kernel to calculate kSq on device
      _computeKSq<D><<<nBlocks, nThreads>>>
            (kSq_.cArray(), minImages_.cArray(), kBasis.cArray(), 
             unitCell().nParameter(), kSize_);
      
      hasKSq_ = true;
   }

   template <int D>
   void WaveList<D>::computedKSq()
   {
      if (hasdKSq_) return; // dKSq already calculated

      // Compute minimum images if needed
      if (!hasMinimumImages_) {
         computeMinimumImages(); 
      }

      // Precondition
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());

      // Calculate dkkBasis and store on device
      int idx;
      HostDArray<cudaReal> dkkBasis_h(unitCell().nParameter() * D * D);
      DeviceArray<cudaReal> dkkBasis;
      for(int i = 0 ; i < unitCell().nParameter(); ++i) {
         for(int j = 0; j < D; ++j) {
            for(int k = 0; k < D; ++k) {
               idx = k + (j * D) + (i * D * D);
               dkkBasis_h[idx] = unitCell().dkkBasis(i, j, k);
            }
         }
      }
      dkkBasis = dkkBasis_h;

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(kSize_, nBlocks, nThreads);

      // Kernel requires block size to be >= the size of dkkBasis.
      // Max size of dkkBasis is 54, so this should always be satisfied
      UTIL_CHECK(nThreads > dkkBasis.capacity()); 

      // Launch kernel to calculate dKSq on device
      size_t sz = sizeof(cudaReal)*dkkBasis.capacity();
      _computedKSq<D><<<nBlocks, nThreads, sz>>>
         (dKSq_.cArray(), minImages_.cArray(), dkkBasis.cArray(), 
          implicitInverse_.cArray(), unitCell().nParameter(), kSize_);
      
      hasdKSq_ = true;
   }

}
}
#endif
