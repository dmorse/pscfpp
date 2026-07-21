#ifndef PRDC_CPU_WAVE_LIST_TPP
#define PRDC_CPU_WAVE_LIST_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"

#include <prdc/field/cpu/FFT.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/hasVariableAngle.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/Sort.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   /*
   * Constructor.
   */
   template <int D>
   WaveList<D>::WaveList(bool isRealField)
    : kSize_(0),
      nBunch_(0),
      isAllocated_(false),
      hasMinImages_(false),
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
   WaveList<D>::~WaveList() 
   {}

   /*
   * Allocate memory used by WaveList.
   */
   template <int D>
   void WaveList<D>::allocate(Mesh<D> const & m, UnitCell<D> const & c) 
   {
      UTIL_CHECK(m.size() > 0);
      UTIL_CHECK(c.nParameter() > 0);
      UTIL_CHECK(!isAllocated_);

      // Create permanent associations with mesh and unit cell
      unitCellPtr_ = &c;
      meshPtr_ = &m;

      // Local copies of properties
      int nParams = unitCell().nParameter();
      IntVec<D> const & meshDimensions = mesh().dimensions();

      // Compute kMeshDimensions_ and kSize_
      if (isRealField_) {
         FFT<D>::computeKMesh(meshDimensions, kMeshDimensions_, kSize_);
      } else {
         kMeshDimensions_ = meshDimensions;
         kSize_ = mesh().size();
      }

      // Allocate memory
      minImages_.allocate(kSize_);
      kSq_.allocate(kMeshDimensions_);
      dKSq_.allocate(nParams);
      for (int i = 0; i < nParams; i++) {
         dKSq_[i].allocate(kMeshDimensions_);
      }

      // Allocate and set up implicitInverse_ array if isRealField_ == true
      // (only depends on mesh dimensions, only used for real fields)
      if (isRealField_) {
         implicitInverse_.allocate(kSize_);

         MeshIterator<D> kItr(kMeshDimensions_);
         int rank;

         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            implicitInverse_[rank] = 
               FFT<D>::hasImplicitInverse(kItr.position(), meshDimensions);
         }
      }

      clearUnitCellData();
      isAllocated_ = true;
   }

   /*
   * Clear all data that depends on unit cell parameters.
   */
   template <int D>
   void WaveList<D>::clearUnitCellData()
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
      }
   }

   /*
   * Compute minimum image vectors for all wavevectors.
   */
   template <int D>
   void WaveList<D>::computeMinimumImages() 
   {
      if (hasMinImages_) return; // min images already calculated

      // Precondition
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());
      UTIL_CHECK(minImages_.capacity() == kSize_);

      MeshIterator<D> kItr(kMeshDimensions_);
      int rank;
      for (kItr.begin(); !kItr.atEnd(); ++kItr) {
         rank = kItr.rank();
         minImages_[rank] = shiftToMinimum(kItr.position(), 
                                           mesh().dimensions(), unitCell());
         kSq_[rank] = unitCell().ksq(minImages_[rank]);
      }

      hasMinImages_ = true;
      hasKSq_ = true;
   }

   /*
   * Compute array of value of |k|^2
   */
   template <int D>
   void WaveList<D>::computeKSq() 
   {
      // If kSq_ is valid, return immediately without recomputing
      if (hasKSq_) return; 

      // If necessary, compute minimum images.
      if (!hasMinImages_) {
         computeMinimumImages(); // computes both min images and kSq
         return;
      }

      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());

      // Compute kSq_
      MeshIterator<D> kItr(kMeshDimensions_);
      int rank;
      for (kItr.begin(); !kItr.atEnd(); ++kItr) {
         rank = kItr.rank();
         kSq_[rank] = unitCell().ksq(minImages_[rank]);
      }
      
      hasKSq_ = true;
   }

   /*
   * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
   */
   template <int D>
   void WaveList<D>::computedKSq()
   {
      if (hasdKSq_) return; // dKSq already calculated

      // Compute minimum images if needed
      if (!hasMinImages_) {
         computeMinimumImages(); 
      }

      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(unitCell().nParameter() > 0);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().isInitialized());

      MeshIterator<D> kItr(kMeshDimensions_);
      int i, rank;
      for (i = 0 ; i < unitCell().nParameter(); ++i) {
         RField<D, CppTp<D> >& dksq = dKSq_[i];
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            dksq[rank] = unitCell().dksq(minImages_[rank], i);
            if (isRealField_) {
               if (implicitInverse_[rank]) {
                  dksq[rank] *= 2.0;
               }
            }
         }
      }
      
      hasdKSq_ = true;
   }

   /*
   * Sort waves by magnitude.
   */
   template <int D>
   void WaveList<D>::sortWaves()
   {
      if (isSorted_) return;

      // Compute wavenumbers if needed
      UTIL_CHECK(isAllocated_);
      if (!hasKSq_) {
         computeKSq();
      }

      // Construct sorted array of items 
      std::vector< Sort::Item<double> > items;
      Sort::Item<double> item;
      for (int i = 0; i < kSize_; ++i) {
         item.value = kSq_[i];
         item.id = i;
         items.push_back(item);
      }
      UTIL_CHECK((int)items.size() == kSize_);
      Sort::sort(items);

      // Fill sortedIds_ array with ids
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

      // Fill bunchIds_ array
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

}
}
}
#endif
