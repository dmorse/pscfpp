#ifndef PSCF_WAVE_LIST_CU_TPP
#define PSCF_WAVE_LIST_CU_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"

#include <pscf/field/cuda/resources.h>
#include <pscf/field/cuda/FFT.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/hasVariableAngle.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/Sort.h>
#include <pscf/math/IntVec.h>

namespace Pscf {

   // CUDA kernels:
   // (defined in an anonymous namespace, used only in this file)
   namespace {

      /*
      * Voronoi-relevant vectors of the lattice of image shifts (see
      * imageShiftNeighbors in shiftToMinimum.h), passed to the kernels
      * by value. Vector j has components v[D*j], ..., v[D*j+D-1].
      */
      struct ImageShifts
      {
         int n;
         int v[14*3];
      };

      /*
      * Square magnitude of the wavevector with integer coordinates img.
      */
      template <int D>
      __device__ cudaReal _kSqImage(int const * img,
                                    cudaReal const * kBasis)
      {
         cudaReal kSq = 0.0;
         for (int c = 0; c < D; ++c) {
            cudaReal kVec = 0.0;
            for (int k = 0; k < D; ++k) {
               kVec += kBasis[D*k + c] * img[k];
            }
            kSq += kVec*kVec;
         }
         return kSq;
      }

      /*
      * Square magnitude of the minimum image r, as stored in kSq by the
      * minimum image kernel (written as in earlier versions of this
      * kernel, so that the stored values keep their rounding).
      */
      template <int D>
      __device__ cudaReal _kSqStored(int const * r,
                                     cudaReal const * kBasis);

      template <>
      __device__ cudaReal _kSqStored<2>(int const * r,
                                        cudaReal const * kBasis)
      {
         cudaReal kVec0 = (kBasis[0] * r[0]) + (kBasis[2] * r[1]);
         cudaReal kVec1 = (kBasis[1] * r[0]) + (kBasis[3] * r[1]);
         return (kVec0*kVec0) + (kVec1*kVec1);
      }

      template <>
      __device__ cudaReal _kSqStored<3>(int const * r,
                                        cudaReal const * kBasis)
      {
         cudaReal kVec0(0.0), kVec1(0.0), kVec2(0.0);
         for (int k = 0; k < 3; ++k) {
            kVec0 += kBasis[3*k]   * r[k];
            kVec1 += kBasis[3*k+1] * r[k];
            kVec2 += kBasis[3*k+2] * r[k];
         }
         return (kVec0*kVec0) + (kVec1*kVec1) + (kVec2*kVec2);
      }

      /*
      * Compute minimum images of each wavevector (D = 2 or 3).
      *
      * On entry, minImages contains the k-grid indices of each
      * wavevector. Each thread computes the minimum images of a set of
      * wavevectors, with the steps of the function minimumImage in
      * shiftToMinimum.cpp: it starts from the image with components in
      * (-meshDims[i]/2, meshDims[i]/2], moves it by the image shift in
      * shifts that lowers kSq most, while one lowers it by more than a
      * fraction epsilon, and then to the largest neighboring image of
      * equal kSq (within a fraction epsilon), repeatedly.
      *
      * When launched, this kernel requires no dynamic memory allocation.
      */
      template <int D> __global__
      void _computeMinimumImages(int* minImages,
                                 cudaReal* kSq,
                                 cudaReal const * kBasis,
                                 int const * meshDims,
                                 int const kSize,
                                 ImageShifts const shifts)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Set epsilon, the relative tolerance for equal kSq (relative,
         // so that the choice of image does not depend on length units)
         #ifdef SINGLE_PRECISION
         cudaReal epsilon = 1.0E-4;
         #else
         cudaReal epsilon = 1.0E-12;
         #endif

         int r[D], u[D];
         for (int i = startID; i < kSize; i += nThreads) {

            // Image with components in (-meshDims[j]/2, meshDims[j]/2]
            for (int j = 0; j < D; ++j) {
               r[j] = minImages[i + j*kSize];
               if (2*r[j] > meshDims[j]) r[j] -= meshDims[j];
            }
            cudaReal kSqMin = _kSqImage<D>(r, kBasis);

            // Move r by the image shift that lowers kSq most
            for (;;) {
               int jBest = -1;
               cudaReal kSqBest = kSqMin;
               for (int j = 0; j < shifts.n; ++j) {
                  for (int k = 0; k < D; ++k) {
                     u[k] = r[k] + shifts.v[D*j + k];
                  }
                  cudaReal g = _kSqImage<D>(u, kBasis);
                  if (g < kSqBest - epsilon*(g > kSqBest ? g : kSqBest)) {
                     jBest = j;
                     kSqBest = g;
                  }
               }
               if (jBest < 0) break;
               for (int k = 0; k < D; ++k) {
                  r[k] += shifts.v[D*jBest + k];
               }
               kSqMin = kSqBest;
            }

            // Move r to a larger (in the order of IntVec) image of
            // equal kSq, for as long as there is one
            bool changed = true;
            while (changed) {
               changed = false;
               for (int j = 0; j < shifts.n; ++j) {
                  for (int k = 0; k < D; ++k) {
                     u[k] = r[k] + shifts.v[D*j + k];
                  }
                  cudaReal g = _kSqImage<D>(u, kBasis);
                  if (fabs(g - kSqMin) <= epsilon*(g > kSqMin ? g : kSqMin)) {
                     bool less = false;  // Is r < u?
                     for (int k = 0; k < D; ++k) {
                        if (r[k] > u[k]) break;
                        if (r[k] < u[k]) {
                           less = true;
                           break;
                        }
                     }
                     if (less) {
                        for (int k = 0; k < D; ++k) {
                           r[k] = u[k];
                        }
                        changed = true;
                     }
                  }
               }
            }

            for (int j = 0; j < D; ++j) {
               minImages[i + j*kSize] = r[j];
            }
            kSq[i] = _kSqStored<D>(r, kBasis);
         }
      }

      /*
      * Compute minimum images of each wavevector, 1D specialization.
      *
      * Kernel must be launched with kSize threads. (One thread calculates
      * one minimum image.) The image shifts are not used in 1D.
      *
      * When launched, this kernel requires no dynamic memory allocation.
      */
      template <> __global__
      void _computeMinimumImages<1>(int* minImages,
                                    cudaReal* kSq,
                                    cudaReal const * kBasis,
                                    int const * meshDims,
                                    int const kSize,
                                    ImageShifts const shifts)
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
      * Compute the kSq array on the GPU.
      */
      template <int D> __global__ 
      void _computeKSq(cudaReal* kSq, 
                       int const * waveBz,
                       cudaReal const * kBasis,
                       int const nParams, 
                       int const kSize)
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
      * Compute the dKSq array on the GPU.
      *
      * The implicitInverse array is only used if isRealField is true.
      * If isRealField is false, implicitInverse may be a null pointer.
      */
      template <int D> __global__ 
      void _computedKSq(cudaReal* dKSq, 
                        int const * waveBz,
                        cudaReal const * dkkBasis,
                        bool const * implicitInverse,
                        int const nParams, 
                        int const kSize,
                        bool isRealField)
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

               if (isRealField) {
                  if (implicitInverse[i]) {
                     // if element i's inverse is implicit
                     dKSqVal *= 2;
                  }
               }

               dKSq[(param * kSize) + i] = dKSqVal;

            } // kSize
         } // nParams
      }

      /*
      * Average dKSq over tied minimum images on the GPU.
      *
      * For each wavevector, finds the images tied with the minimum image
      * in waveBz as tiedMinimumImages in the CPU WaveList does: the
      * images reached by steps of the vectors in shifts through images
      * with kSq equal to that of the minimum image within a fraction
      * epsilon. If there are several, dKSq is set as by tiedImageDKSq
      * to the average of their derivatives, unless these all agree
      * within this fraction, and multiplied by 2 if the wave has an
      * implicit inverse (if isRealField is true). If a wave has more
      * than 2^D tied images, *overflow is set to 1.
      *
      * When launched, this kernel requires dynamic memory allocation of
      * sizeof(cudaReal) * nParams * D * D.
      */
      template <int D> __global__
      void _averageTiedDKSq(cudaReal* dKSq,
                            int const * waveBz,
                            cudaReal const * dkkBasis,
                            cudaReal const * kBasis,
                            ImageShifts const shifts,
                            bool const * implicitInverse,
                            int const nParams,
                            int const kSize,
                            bool isRealField,
                            int* overflow)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;

         // Load dkkBasis into shared memory for fast access
         extern __shared__ cudaReal dkkBasis_[];
         if (threadIdx.x < nParams * D * D) {
            dkkBasis_[threadIdx.x] = dkkBasis[threadIdx.x];
         }
         __syncthreads(); // wait for all threads to finish

         // Relative tolerance, as in the CPU WaveList
         #ifdef SINGLE_PRECISION
         cudaReal epsilon = 1.0E-4;
         #else
         cudaReal epsilon = 1.0E-9;
         #endif

         // Tied images, image t with components tied[D*t], ... (a
         // minimum image is tied with at most 2^D images)
         const int maxTied = 1 << D;
         int tied[maxTied * D];
         int u[D];

         for (int i = startID; i < kSize; i += nThreads) {

            // Find the images tied with the minimum image
            for (int j = 0; j < D; ++j) {
               tied[j] = waveBz[i + (j * kSize)];
            }
            int nTied = 1;
            cudaReal g0 = _kSqImage<D>(tied, kBasis);
            for (int t = 0; t < nTied; ++t) {
               for (int s = 0; s < shifts.n; ++s) {
                  for (int j = 0; j < D; ++j) {
                     u[j] = tied[D*t + j] + shifts.v[D*s + j];
                  }
                  cudaReal g = _kSqImage<D>(u, kBasis);
                  if (fabs(g - g0) > epsilon*(g > g0 ? g : g0)) continue;
                  bool found = false;
                  for (int m = 0; m < nTied; ++m) {
                     bool same = true;
                     for (int j = 0; j < D; ++j) {
                        if (tied[D*m + j] != u[j]) same = false;
                     }
                     if (same) found = true;
                  }
                  if (found) continue;
                  if (nTied == maxTied) {
                     *overflow = 1;
                     continue;
                  }
                  for (int j = 0; j < D; ++j) {
                     tied[D*nTied + j] = u[j];
                  }
                  ++nTied;
               }
            }
            if (nTied == 1) continue;

            // Average the derivatives over the tied images, unless they
            // all agree (then use that of the first, as _computedKSq)
            for (int param = 0; param < nParams; ++param) {
               cudaReal value0 = 0.0;
               cudaReal sum = 0.0;
               cudaReal maxDiff = 0.0;
               cudaReal scale = 0.0;
               for (int t = 0; t < nTied; ++t) {
                  cudaReal value = 0.0;
                  for (int j = 0; j < D; ++j) {
                     for (int k = 0; k < D; ++k) {
                        value += cudaReal(tied[D*t + j])
                                 * cudaReal(tied[D*t + k])
                                 * dkkBasis_[k + (j * D) + (param * D * D)];
                     }
                  }
                  if (t == 0) value0 = value;
                  sum += value;
                  if (fabs(value - value0) > maxDiff) {
                     maxDiff = fabs(value - value0);
                  }
                  if (fabs(value) > scale) scale = fabs(value);
               }
               cudaReal dKSqVal = value0;
               if (maxDiff > epsilon*scale) {
                  dKSqVal = sum/cudaReal(nTied);
               }
               if (isRealField) {
                  if (implicitInverse[i]) {
                     dKSqVal *= 2;
                  }
               }
               dKSq[(param * kSize) + i] = dKSqVal;
            }

         } // kSize
      }

   } // anonymous namespace Pscf::Prdc::(unnamed)

   // Member function definitions

   /*
   * Constructor.
   */
   template <int D>
   WaveList<D,CUT>::WaveList(bool isRealField)
    : kSize_(0),
      nBunch_(0),
      isAllocated_(false),
      hasMinImages_(false),
      hasMinImages_h_(false),
      hasKSq_(false),
      hasdKSq_(false),
      isSorted_(false),
      isRealField_(isRealField),
      unitCellPtr_(nullptr),
      meshPtr_(nullptr)
   {}

   /*
   * Destructor.
   */
   template <int D>
   WaveList<D,CUT>::~WaveList()
   {

      if (dKSqSlices_.isAllocated()) {
         int nParams = dKSqSlices_.capacity();
         for (int i = 0; i < nParams; i++) {
            if (dKSqSlices_[i].isAssociated()) {
               dKSqSlices_[i].dissociate();
            }
         }
      }

      /*
      * The above loop dissociates elements of dKSqSlices_ from memory
      * owned by dKSq_. Because the destructor body is called before
      * destructors for members, this operation will occur before either
      * container is destroyed.  Strictly speaking, this is a redundant
      * safety measure as long as dKSqSlices_ appears in the declaration
      * of members after dKSq_, because this order causes dKSqSlices_
      * to to be destroyed before dKSq_, and the destructor for each
      * element of dKSqSlices_ would also call dissociate() if the
      * association still existed at that point.
      */

   }

   /*
   * Allocate memory.
   *
   * This function allocates and construct implicitInverse_d_ if and
   * only if isRealField_ == true.
   */
   template <int D>
   void WaveList<D,CUT>::allocate(Mesh<D> const & m, UnitCell<D> const & c)
   {
      UTIL_CHECK(m.size() > 0);
      UTIL_CHECK(c.nParameter() > 0);
      UTIL_CHECK(!isAllocated_);

      // Create permanent associations with mesh and unit cell
      unitCellPtr_ = &c;
      meshPtr_ = &m;

      int nParams = unitCell().nParameter();
      IntVec<D> const & meshDimensions = mesh().dimensions();

      // Compute kMeshDimensions_ and kSize_
      if (isRealField_) {
         FFT<D,CUT>::computeKMesh(meshDimensions, kMeshDimensions_, kSize_);
      } else {
         kMeshDimensions_ = meshDimensions;
         kSize_ = mesh().size();
      }
      UTIL_CHECK(kSize_ > 0);

      minImages_.allocate(kSize_ * D);
      kSq_.allocate(kMeshDimensions_);
      dKSq_.allocate(kSize_ * nParams);
      dKSqSlices_.allocate(nParams);
      for (int i = 0; i < nParams; i++) {
         dKSqSlices_[i].associate(dKSq_, i*kSize_, kMeshDimensions_);
      }

      // Allocate and fill implicitInverse arrays if isRealField_ == true
      // (only depends on mesh dimensions, only used for real fields)
      if (isRealField_) {
         implicitInverse_d_.allocate(kSize_);
         implicitInverse_h_.allocate(kSize_);
         HostSendArray<bool,CUT> implicitTemp;
         implicitTemp.allocate(kSize_);
         MeshIterator<D> kItr(kMeshDimensions_);
         int rank;
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            implicitTemp[rank] = 
               FFT<D,CUT>::hasImplicitInverse(kItr.position(), 
                                              meshDimensions);
            implicitInverse_h_[rank] = implicitTemp[rank];
         }
         implicitInverse_d_ = implicitTemp; // copy to device memory
      }

      clearUnitCellData();
      isAllocated_ = true;
   }

   /*
   * Clear data that depends on unit cell parameters.
   */
   template <int D>
   void WaveList<D,CUT>::clearUnitCellData()
   {
      hasKSq_ = false;
      hasdKSq_ = false;
      if (unitCell().nParameter() > 1) {
         isSorted_ = false;
         sortedBunches_.clear();
         nBunch_ = 0;
      }
      if (hasVariableAngle<D>(unitCell().lattice())) {
         hasMinImages_ = false;
         hasMinImages_h_ = false;
      }
   }

   /*
   * Compute minimum image vectors and kSq.
   */
   template <int D>
   void WaveList<D,CUT>::computeMinimumImages()
   {
      // If min images are valid, return immediately
      if (hasMinImages_) return;

      // Precondition
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(minImages_.capacity() == kSize_ * D);

      // Set initial array of images to contain the k-grid points
      HostSendArray<int,CUT> imagesTmp(D*kSize_);
      MeshIterator<D> kItr(kMeshDimensions_);
      for (int i = 0; i < D; i++) {
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            imagesTmp[kItr.rank() + (i*kSize_)] = kItr.position(i);
         }
      }
      minImages_ = imagesTmp; // copy to device

      // Get kBasis and meshDims and store on device
      HostSendArray<cudaReal,CUT> kBasis_h(D*D);
      HostSendArray<int,CUT> meshDims_h(D);
      int idx = 0;
      for (int j = 0; j < D; ++j) {
         for (int k = 0; k < D; ++k) {
            kBasis_h[idx] = unitCell().kBasis(j)[k];
            idx++;
         }
         meshDims_h[j] = mesh().dimension(j);
      }

      BackendArray<cudaReal,CUT> kBasis(D*D);
      BackendArray<int,CUT> meshDims(D);
      kBasis = kBasis_h;
      meshDims = meshDims_h;

      // Voronoi-relevant vectors of the lattice of image shifts
      FSArray<IntVec<D>, 14> neighbors;
      imageShiftNeighbors<D>(mesh().dimensions(), unitCell(), neighbors);
      ImageShifts shifts;
      shifts.n = neighbors.size();
      for (int j = 0; j < neighbors.size(); ++j) {
         for (int i = 0; i < D; ++i) {
            shifts.v[D*j + i] = neighbors[j][i];
         }
      }

      // GPU resources (one thread per wavevector)
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(kSize_, nBlocks, nThreads);

      if ((D == 3) && (nThreads < 128)) {
         // Thread blocks too small. Manually set nThreads to 128
         // (computedKSq needs more threads per block than dkkBasis has
         // elements, up to 54 in 3D)
         ThreadArray::setThreadsPerBlock(128);
         ThreadArray::setThreadsLogical(kSize_, nBlocks, nThreads);

         // If the above was successful, print warning
         Log::file() << "\nWarning: "
                     << "nThreads too small for computeMinimumImages.\n"
                     << "Setting nThreads equal to 128." << std::endl;
      }

      UTIL_CHECK(kSq_.isAllocated());
      UTIL_CHECK(minImages_.isAllocated());

      // Launch kernel
      _computeMinimumImages<D><<<nBlocks, nThreads>>>
         (minImages_.cArray(), kSq_.cArray(), kBasis.cArray(),
          meshDims.cArray(), kSize_, shifts);

      hasMinImages_ = true;
      hasKSq_ = true;
   }

   /*
   * Compute values of k^2, using existing minImages if possible.
   */
   template <int D>
   void WaveList<D,CUT>::computeKSq()
   {
      // If kSq values are valid, return immediately without recomputing
      if (hasKSq_) return;

      // If necessary, compute minimum images
      if (!hasMinImages_) {
         computeMinimumImages(); // compute both min images and kSq
         return;
      }

      // If this point is reached, calculate kSq using _computeKSq kernel

      // Preconditions
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(isAllocated_);

      // Get kBasis and store on device
      HostSendArray<cudaReal,CUT> kBasis_h(D*D);
      BackendArray<cudaReal,CUT> kBasis(D*D);
      int idx = 0;
      for (int j = 0; j < D; ++j) {
         for (int k = 0; k < D; ++k) {
            kBasis_h[idx] = unitCell().kBasis(j)[k];
            idx++;
         }
      }
      kBasis = kBasis_h;

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(kSize_, nBlocks, nThreads);

      // Launch kernel to calculate kSq on device
      UTIL_CHECK(kSq_.isAllocated());
      UTIL_CHECK(minImages_.isAllocated());
      _computeKSq<D><<<nBlocks, nThreads>>>
            (kSq_.cArray(), minImages_.cArray(), kBasis.cArray(),
             unitCell().nParameter(), kSize_);

      hasKSq_ = true;
   }

   /*
   * Compute derivatives of k^2 with respect to unit cell parameters.
   */
   template <int D>
   void WaveList<D,CUT>::computedKSq()
   {
      if (hasdKSq_) return; // dKSq already calculated

      // Compute minimum images if needed
      if (!hasMinImages_) {
         computeMinimumImages();
      }

      // Preconditions
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(isAllocated_);

      // Calculate dkkBasis and store on device
      int idx;
      HostSendArray<cudaReal,CUT> dkkBasis_h(unitCell().nParameter() * D * D);
      BackendArray<cudaReal,CUT> dkkBasis;
      for (int i = 0 ; i < unitCell().nParameter(); ++i) {
         for (int j = 0; j < D; ++j) {
            for (int k = 0; k < D; ++k) {
               idx = k + (j * D) + (i * D * D);
               dkkBasis_h[idx] = unitCell().dkkBasis(i, j, k);
            }
         }
      }
      dkkBasis = dkkBasis_h;

      // GPU resources
      int nBlocks, nThreads;
      ThreadArray::setThreadsLogical(kSize_, nBlocks, nThreads);

      // Kernel requires block size to be >= the size of dkkBasis.
      // Max size of dkkBasis is 54, so this should always be satisfied
      UTIL_CHECK(nThreads > dkkBasis.capacity());

      // Preconditions for kernel
      UTIL_CHECK(dKSq_.isAllocated());
      UTIL_CHECK(minImages_.isAllocated());
      bool const * implicitInversePtr = nullptr;
      if (isRealField_) {
         UTIL_CHECK(implicitInverse_d_.isAllocated());
         implicitInversePtr = implicitInverse_d_.cArray();
      }

      // Launch kernel to calculate dKSq on device
      size_t sz = sizeof(cudaReal)*dkkBasis.capacity();
      _computedKSq<D><<<nBlocks, nThreads, sz>>>
         (dKSq_.cArray(), minImages_.cArray(), dkkBasis.cArray(),
          implicitInversePtr, unitCell().nParameter(), kSize_,
          isRealField_);

      // Minimum images tied with images whose |k|^2 has other derivatives
      // (cells with a variable angle only, e.g. rhombic at 120 degrees):
      // use the average over the tied images, as the CPU WaveList does
      if (hasVariableAngle<D>(unitCell().lattice())) {
         HostSendArray<cudaReal,CUT> kBasis_h(D*D);
         BackendArray<cudaReal,CUT> kBasis(D*D);
         for (int j = 0; j < D; ++j) {
            for (int k = 0; k < D; ++k) {
               kBasis_h[D*j + k] = unitCell().kBasis(j)[k];
            }
         }
         kBasis = kBasis_h;
         FSArray<IntVec<D>, 14> neighbors;
         imageShiftNeighbors<D>(mesh().dimensions(), unitCell(), neighbors);
         ImageShifts shifts;
         shifts.n = neighbors.size();
         for (int j = 0; j < neighbors.size(); ++j) {
            for (int k = 0; k < D; ++k) {
               shifts.v[D*j + k] = neighbors[j][k];
            }
         }
         HostSendArray<int,CUT> overflow_h(1);
         overflow_h[0] = 0;
         BackendArray<int,CUT> overflow(1);
         overflow = overflow_h;
         _averageTiedDKSq<D><<<nBlocks, nThreads, sz>>>
            (dKSq_.cArray(), minImages_.cArray(), dkkBasis.cArray(),
             kBasis.cArray(), shifts, implicitInversePtr,
             unitCell().nParameter(), kSize_, isRealField_,
             overflow.cArray());
         HostRecvArray<int,CUT> overflowOut;
         overflowOut = overflow;
         if (overflowOut[0]) {
            UTIL_THROW("More than 2^D tied minimum images");
         }
      }

      hasdKSq_ = true;
   }

   /*
   * Sort waves by magnitude.
   */
   template <int D>
   void WaveList<D,CUT>::sortWaves()
   {
      // If waves are already sorted, return immediately
      if (isSorted_) return;

      // Compute wavenumbers if necessary
      UTIL_CHECK(isAllocated_);
      if (!hasKSq_) {
         computeKSq();
      }

      // Copy values of kSq to host
      HostRecvArray<cudaReal,CUT> kSq_h(kSq_);

      // Construct Sort::Item objects with value = kSq, id = wave id
      std::vector< Sort::Item<double> > items;
      Sort::Item<double> item;
      for (int i = 0; i < kSize_; ++i) {
         item.value = kSq_h[i];
         item.id = i;
         items.push_back(item);
      }
      UTIL_CHECK((int)items.size() == kSize_);

      // Sort array of items
      Sort::sort(items);

      // Fill sortedIds_ member array with sorted ids
      if (!sortedIds_.isAllocated()) {
         sortedIds_.allocate(kSize_);
      }
      UTIL_CHECK(sortedIds_.capacity() == kSize_);
      for (int i = 0; i < kSize_; ++i) {
         sortedIds_[i] = items[i].id;
      }

      // Construct sortedBunches_ array and set nBunch_
      double epsilon = 1.0E-8;
      sortedBunches_.clear();
      Sort::findBunches(items, sortedBunches_, epsilon);
      nBunch_ = sortedBunches_.size();

      // Construct bunchIds_ array
      UTIL_CHECK(nBunch_ > 0);
      if (!bunchIds_.isAllocated()) {
         bunchIds_.allocate(kSize_);
      }
      UTIL_CHECK(bunchIds_.capacity() == kSize_);
      int begin, end, ib, iw;
      for (ib = 0; ib < nBunch_; ++ib) {
         begin = sortedBunches_[ib][0];
         end = sortedBunches_[ib][1];
         UTIL_CHECK(end > begin);
         for (iw = begin; iw < end; ++iw) {
            bunchIds_[sortedIds_[iw]] = ib;
         }
      }
      UTIL_CHECK(end == kSize_);

      isSorted_ = true;
   }

   /*
   * Get minimum images on the host by reference.
   */
   template <int D>
   DArray< IntVec<D> > const & 
   WaveList<D,CUT>::minImages_h() const
   {
      UTIL_CHECK(hasMinImages_);
      if (!hasMinImages_h_) {
         if (!minImages_h_.isAllocated()) {
            minImages_h_.allocate(kSize_);
         }
         HostRecvArray<int,CUT> minImages_temp;
         minImages_temp = minImages_;
         int offset;
         for (int j = 0; j < D; ++j) {
            offset = j*kSize_;
            for (int i = 0; i < kSize_; ++i) {
               minImages_h_[i][j] = minImages_temp[i + offset];
            }
         }
         hasMinImages_h_ = true;
      }
      return minImages_h_;
   }

} // namespace Pscf
#endif
