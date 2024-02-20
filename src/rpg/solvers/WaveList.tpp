#ifndef PSPG_WAVE_LIST_TPP
#define PSPG_WAVE_LIST_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"
#include "cuComplex.h"
#include <pscf/cuda/GpuResources.h>

namespace Pscf {
namespace Rpg
{

   // Need a reference table that maps index to a pair wavevector
   // Ideally we can have a group of thread dealing with only
   // the non-implicit part and the implicit part
   static __global__ 
   void makeDksqHelperWave(cudaReal* dksq, const int* waveBz,
                           const cudaReal* dkkBasis,
                           const int* partnerId,
                           const int* selfId,
                           const bool* implicit,
                           int nParams, int kSize,
                                             int size, int dim) 
   {
      // Actual size is nStar*nParams
      // Each thread does nParams calculation
      // Big performance hit if thread >= 0.5dimension(n-1)
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      
      //loop through the entire array
      int pId;
      for (int param = 0; param < nParams; ++param) {
         for (int i = startID; i < size; i += nThreads) {
            for (int j = 0; j < dim; ++j) {
               for (int k = 0; k < dim; ++k) {
                  if (!implicit[i]) {
                     // Not = need += so still need memset
                     dksq[(param * size) + i] 
                         += waveBz[selfId[i] * dim + j] 
                            * waveBz[ selfId[i] * dim + k]
                            * dkkBasis[k + (j * dim) + (param * dim * dim)];
                  } else {
                     pId = partnerId[i];
                     dksq[(param * size) + i] += waveBz[selfId[pId] * dim + j]
                        * waveBz[selfId[pId] * dim + k]
                        * dkkBasis[k + (j * dim) + (param * dim * dim)];
                  }
               } //dim
            } //dim
         } //size
      } //nParams
   }

   static __global__ 
   void makeDksqReduction(cudaReal* dksq, const int* partnerId,
                          int nParams, int kSize, int rSize) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;

      // Add i in the implicit part into their partner's result
      int pId;
      for(int param = 0; param < nParams; ++param) {
         for (int i = startID + kSize; i < rSize; i += nThreads) {
            pId = partnerId[i];
            dksq[(param * rSize) + pId] += dksq[(param * rSize) + i];
         }
      }
   }

   template <int D>
   WaveList<D>::WaveList()
    : minImage_d(nullptr),
      dkSq_(nullptr),
      partnerIdTable(nullptr),
      partnerIdTable_d(nullptr),
      kSize_(0),
      rSize_(0),
      nParams_(0),
      isAllocated_(false),
      hasMinimumImages_(false)
   {}

   template <int D>
   WaveList<D>::~WaveList() {
      if (isAllocated_) {
         cudaFree(minImage_d);
         cudaFree(dkSq_);
         cudaFree(partnerIdTable_d);
         cudaFree(selfIdTable_d);
         cudaFree(implicit_d);
         cudaFree(dkkBasis_d);
         delete[] kSq_;
         delete implicit;
      }
   }

   template <int D>
   void WaveList<D>::allocate(Mesh<D> const & mesh, 
                              UnitCell<D> const & unitCell) 
   {
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(unitCell.nParameter() > 0);
      UTIL_CHECK(!isAllocated_);

      rSize_ = mesh.size();
      dimensions_ = mesh.dimensions();
      nParams_ = unitCell.nParameter();

      // Compute DFT mesh size kSize_
      kSize_ = 1;
      for(int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kSize_ *= mesh.dimension(i);
         } else {
            kSize_ *= (mesh.dimension(i)/ 2 + 1);
         }
      }

      minImage_.allocate(kSize_);
      gpuErrchk(
         cudaMalloc((void**) &minImage_d, sizeof(int) * rSize_ * D)
      );

      kSq_ = new cudaReal[rSize_];
      gpuErrchk(
         cudaMalloc((void**) &dkSq_, sizeof(cudaReal) * rSize_ * nParams_)
      );

      partnerIdTable = new int[mesh.size()];
      gpuErrchk(
         cudaMalloc((void**) &partnerIdTable_d, sizeof(int) * mesh.size())
      );

      selfIdTable = new int[mesh.size()];
      gpuErrchk(
         cudaMalloc((void**) &selfIdTable_d, sizeof(int) * mesh.size())
      );

      implicit = new bool[mesh.size()];
      gpuErrchk(
         cudaMalloc((void**) &implicit_d, sizeof(bool) * mesh.size())
      );

      dkkBasis = new cudaReal[6 * D * D];
      gpuErrchk(
         cudaMalloc((void**) &dkkBasis_d, sizeof(cudaReal) * 6 * D * D)
      );

      isAllocated_ = true;
   }

   template <int D>
   void WaveList<D>::computeMinimumImages(Mesh<D> const & mesh, 
                                          UnitCell<D> const & unitCell) {
      // Precondition
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(unitCell.nParameter() > 0);
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.isInitialized());

      MeshIterator<D> itr(mesh.dimensions());
      IntVec<D> waveId;
      IntVec<D> G2;
      IntVec<D> tempIntVec;
      int partnerId;

      //min image needs mesh size of them
      //partner only need kSize of them
      //does setting iterator over kdim solves thing?
      int kDimRank = 0;
      int implicitRank = kSize_;
      //kDimRank + implicitRank = rSize
      int* invertedIdTable = new int[rSize_];

      for (itr.begin(); !itr.atEnd(); ++itr) {         
         // If not implicit
         if (itr.position(D - 1) < mesh.dimension(D-1)/2 + 1) {
            implicit[kDimRank] = false;
            selfIdTable[kDimRank] = itr.rank();
            invertedIdTable[itr.rank()] = kDimRank;
            kDimRank++;
         } else {
            implicit[implicitRank] = true;
            selfIdTable[implicitRank] = itr.rank();
            invertedIdTable[itr.rank()] = implicitRank;
            implicitRank++;
         }
      }

      int* tempMinImage = new int[rSize_ * D];
      for (itr.begin(); !itr.atEnd(); ++itr) {
         kSq_[itr.rank()] = unitCell.ksq(itr.position());

         #if 0
         //we get position but set mesh dim to be larger, should be okay
         shiftToMinimum(itr.position(), mesh.dimensions(), 
                        minImage_ + (itr.rank() * D));
         #endif

         // We get position but set mesh dim to be larger, should be okay
         // not the most elegant code with repeated copying but reduces 
         // repeated code from pscf
         waveId = itr.position();
         tempIntVec = shiftToMinimum(waveId, mesh.dimensions(), unitCell);
         for(int i = 0; i < D; i++) {
            (tempMinImage + (itr.rank() * D))[i] = tempIntVec[i];
         }
         
         if(itr.position(D - 1) < mesh.dimension(D-1)/2 + 1) {
            minImage_[invertedIdTable[itr.rank()]] = tempIntVec;
         }

         for(int j = 0; j < D; ++j) {
            G2[j] = -waveId[j];
         }
         mesh.shift(G2);
         partnerId = mesh.rank(G2);
         partnerIdTable[invertedIdTable[itr.rank()]] = invertedIdTable[partnerId];
      }

      #if 0
      /*
      std::cout<< "Sum kDimRank implicitRank: "
                << kDimRank + (implicitRank-kSize_)<<std::endl;
      std::cout<< "This is kDimRank sanity check "<<kDimRank<<std::endl;
      for (int i = 0; i < rSize_; i++) {
         std::cout << i << ' '
                   << selfIdTable[i]<< ' '<<partnerIdTable[i]<<' '
                   << implicit[i] << std::endl;
         }
      */
      #endif

      gpuErrchk(
         cudaMemcpy(minImage_d, tempMinImage, 
                    sizeof(int) * rSize_ * D, cudaMemcpyHostToDevice)
      );

      // Partner is much smaller but we keep this for now
      gpuErrchk(
         cudaMemcpy(partnerIdTable_d, partnerIdTable, 
                    sizeof(int) * mesh.size(), cudaMemcpyHostToDevice)
      );

      gpuErrchk(
         cudaMemcpy(selfIdTable_d, selfIdTable, 
                    sizeof(int) * mesh.size(), cudaMemcpyHostToDevice)
      );

      gpuErrchk(
         cudaMemcpy(implicit_d, implicit, 
                    sizeof(bool) * mesh.size(), cudaMemcpyHostToDevice)
      );
      
      delete[] tempMinImage;
      hasMinimumImages_ = true;
   }

   template <int D>
   void WaveList<D>::computeKSq(UnitCell<D> const & unitCell) {
      //pass for now
   }

   template <int D>
   void WaveList<D>::computedKSq(UnitCell<D> const & unitCell)
   {
      // dkkbasis is something determined from unit cell size
      // min image needs to be on device but okay since its only done once.
      // Second to last parameter is number of stars originally

      // Precondition
      UTIL_CHECK(hasMinimumImages_);
      UTIL_CHECK(unitCell.nParameter() > 0);
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.isInitialized());

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(rSize_, nBlocks, nThreads);

      int idx;
      for(int i = 0 ; i < unitCell.nParameter(); ++i) {
         for(int j = 0; j < D; ++j) {
            for(int k = 0; k < D; ++k) {
               idx = k + (j * D) + (i * D * D);
               dkkBasis[idx] = unitCell.dkkBasis(i, j, k);
            }
         }
      }

      cudaMemcpy(dkkBasis_d, dkkBasis,
                 sizeof(cudaReal) * unitCell.nParameter() * D * D,
                 cudaMemcpyHostToDevice);
      cudaMemset(dkSq_, 0, 
                 unitCell.nParameter() * rSize_ * sizeof(cudaReal));

      makeDksqHelperWave<<<nBlocks, nThreads>>>
         (dkSq_, minImage_d, dkkBasis_d, partnerIdTable_d,
          selfIdTable_d, implicit_d, unitCell.nParameter(), 
          kSize_, rSize_, D);
       
      makeDksqReduction<<<nBlocks, nThreads>>>
          (dkSq_, partnerIdTable_d, unitCell.nParameter(),
            kSize_, rSize_);
   }


}
}
#endif
