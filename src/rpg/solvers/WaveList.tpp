#ifndef RPG_WAVE_LIST_TPP
#define RPG_WAVE_LIST_TPP

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
                     dksq[(param * size) + i] += 
                        waveBz[selfId[pId] * dim + j]
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
    : kSize_(0),
      rSize_(0),
      nParams_(0),
      isAllocated_(false),
      hasMinimumImages_(false)
   {}

   template <int D>
   WaveList<D>::~WaveList() 
   {}

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

      minImage_h_.allocate(kSize_);
      minImage_.allocate(rSize_ * D);

      kSq_h_.allocate(rSize_);
      dkSq_.allocate(rSize_ * nParams_);

      partnerIdTable_h_.allocate(mesh.size());
      partnerIdTable_.allocate(mesh.size());

      selfIdTable_h_.allocate(mesh.size());
      selfIdTable_.allocate(mesh.size());

      implicit_h_.allocate(mesh.size());
      implicit_.allocate(mesh.size());

      dkkBasis_h_.allocate(nParams_ * D * D);
      dkkBasis_.allocate(nParams_ * D * D);

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
      HostDArray<int> invertedIdTable(rSize_);

      for (itr.begin(); !itr.atEnd(); ++itr) {         
         // If not implicit
         if (itr.position(D - 1) < mesh.dimension(D-1)/2 + 1) {
            implicit_h_[kDimRank] = false;
            selfIdTable_h_[kDimRank] = itr.rank();
            invertedIdTable[itr.rank()] = kDimRank;
            kDimRank++;
         } else {
            implicit_h_[implicitRank] = true;
            selfIdTable_h_[implicitRank] = itr.rank();
            invertedIdTable[itr.rank()] = implicitRank;
            implicitRank++;
         }
      }

      HostDArray<int> tempMinImage(rSize_ * D);
      for (itr.begin(); !itr.atEnd(); ++itr) {
         kSq_h_[itr.rank()] = unitCell.ksq(itr.position());

         // We get position but set mesh dim to be larger, should be okay
         // not the most elegant code with repeated copying but reduces 
         // repeated code from pscf
         waveId = itr.position();
         tempIntVec = shiftToMinimum(waveId, mesh.dimensions(), unitCell);
         for(int i = 0; i < D; i++) {
            tempMinImage[itr.rank() * D + i] = tempIntVec[i];
         }
         
         if(itr.position(D - 1) < mesh.dimension(D-1)/2 + 1) {
            minImage_h_[invertedIdTable[itr.rank()]] = tempIntVec;
         }

         for(int j = 0; j < D; ++j) {
            G2[j] = -waveId[j];
         }
         mesh.shift(G2);
         partnerId = mesh.rank(G2);
         partnerIdTable_h_[invertedIdTable[itr.rank()]] = 
                                             invertedIdTable[partnerId];
      }

      // Transfer to device
      minImage_ = tempMinImage;
      selfIdTable_ = selfIdTable_h_;
      implicit_ = implicit_h_;
      // Partner is much smaller but we keep this for now
      partnerIdTable_ = partnerIdTable_h_;
      
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
               dkkBasis_h_[idx] = unitCell.dkkBasis(i, j, k);
            }
         }
      }

      dkkBasis_ = dkkBasis_h_; // transfer dkkBasis to device

      VecOp::eqS(dkSq_, 0); // initialize dkSq_ to an array of 0s

      makeDksqHelperWave<<<nBlocks, nThreads>>>
         (dkSq_.cArray(), minImage_.cArray(), dkkBasis_.cArray(), 
          partnerIdTable_.cArray(), selfIdTable_.cArray(), 
          implicit_.cArray(), unitCell.nParameter(), kSize_, rSize_, D);
       
      makeDksqReduction<<<nBlocks, nThreads>>>
          (dkSq_.cArray(), partnerIdTable_.cArray(), 
           unitCell.nParameter(), kSize_, rSize_);
   }


}
}
#endif
