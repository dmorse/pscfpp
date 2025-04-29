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

      FFT<D>::computeKMesh(mesh().dimensions(), kMeshDimensions_, kSize_);

      minImages_.allocate(kSize_);
      kSq_.allocate(kMeshDimensions_);
      dKSq_.allocate(nParams);
      for (int i = 0; i < nParams; i++) {
         dKSq_[i].allocate(kMeshDimensions_);
      }
      implicitInverse_.allocate(kSize_);

      // Set up implicitInverse_ array (only depends on mesh dimensions)
      MeshIterator<D> kItr(kMeshDimensions_);
      int inverseId;
      for (kItr.begin(); !kItr.atEnd(); ++kItr) {
         if (kItr.position(D-1) == 0) {
            inverseId = 0;
         } else {
            inverseId = mesh().dimension(D-1) - kItr.position(D-1);
         }
         if (inverseId > kMeshDimensions_[D-1]) {
            implicitInverse_[kItr.rank()] = true;
         } else {
            implicitInverse_[kItr.rank()] = false;
         }
      }

      isAllocated_ = true;
   }

   template <int D>
   void WaveList<D>::clearUnitCellData()
   {
      hasKSq_ = false;
      hasdKSq_ = false;
      if (hasVariableAngle<D>(unitCell().lattice())) {
         hasMinimumImages_ = false;
      }
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

   template <int D>
   void WaveList<D>::computeKSq() 
   {
      if (hasKSq_) return; // kSq already calculated

      if (!hasMinimumImages_) {
         computeMinimumImages(); // compute both min images and kSq
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
         RField<D>& field = dKSq_[i];
         for (kItr.begin(); !kItr.atEnd(); ++kItr) {
            rank = kItr.rank();
            field[rank] = unitCell().dksq(minImages_[rank], i);
         }
      }
      
      hasdKSq_ = true;
   }

}
}
}
#endif
