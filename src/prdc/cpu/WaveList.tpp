#ifndef PRDC_CPU_WAVE_LIST_TPP
#define PRDC_CPU_WAVE_LIST_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"

#include <prdc/cpu/FFT.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/hasVariableAngle.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   /*
   * Constructor.
   */
   template <int D>
   WaveList<D>::WaveList(bool isRealField)
    : kSize_(0),
      isAllocated_(false),
      hasMinimumImages_(false),
      hasKSq_(false),
      hasdKSq_(false),
      unitCellPtr_(nullptr),
      meshPtr_(nullptr)
   {  isRealField_ = isRealField; }

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
      if (hasVariableAngle<D>(unitCell().lattice())) {
         hasMinimumImages_ = false;
      }
   }

   /*
   * Compute minimum image vectors for all wavevectors.
   */
   template <int D>
   void WaveList<D>::computeMinimumImages() 
   {
      if (hasMinimumImages_) return; // min images already calculated

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

      hasMinimumImages_ = true;
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
      if (!hasMinimumImages_) {
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
      if (!hasMinimumImages_) {
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
         RField<D>& dksq = dKSq_[i];
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

}
}
}
#endif
