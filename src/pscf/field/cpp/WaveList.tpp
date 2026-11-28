#ifndef PSCF_WAVE_LIST_CP_TPP
#define PSCF_WAVE_LIST_CP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.h"

#include <pscf/field/cpp/FFT.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/crystal/hasVariableAngle.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/Sort.h>

#include <algorithm>
#include <cmath>

namespace Pscf {

   // Functions used to average dKSq over tied minimum images
   // (defined in an anonymous namespace, used only in this file)
   namespace {

      /*
      * Minimum images tied with the minimum image r.
      *
      * A minimum image on the boundary of the Voronoi cell of the
      * lattice of image shifts is tied with other images of the same
      * wavevector. This function returns in tied the image r (first)
      * and the images with |k|^2 equal to that of r within a fraction
      * 1.0E-9 that are reached from r by steps of the vectors in
      * neighbors (from imageShiftNeighbors) through such images. These
      * steps connect all tied minimum images, of which there are at
      * most 2^D. The tolerance is larger than that of shiftToMinimum,
      * so that cell parameters written to a file (with 11 significant
      * digits) at a geometry with tied images give the same images.
      */
      template <int D>
      void tiedMinimumImages(IntVec<D> const & r,
                             FSArray<IntVec<D>, 14> const & neighbors,
                             UnitCell<D> const & cell,
                             FSArray<IntVec<D>, 8>& tied)
      {
         const double epsilon = 1.0E-9;
         const double g0 = cell.ksq(r);
         tied.clear();
         tied.append(r);
         IntVec<D> u;
         for (int t = 0; t < tied.size(); ++t) {
            for (int j = 0; j < neighbors.size(); ++j) {
               u.add(tied[t], neighbors[j]);
               double g = cell.ksq(u);
               if (std::abs(g - g0) > epsilon*std::max(g, g0)) continue;
               bool found = false;
               for (int m = 0; m < tied.size(); ++m) {
                  if (tied[m] == u) found = true;
               }
               if (!found) {
                  if (tied.size() == (1 << D)) {
                     UTIL_THROW("More than 2^D tied minimum images");
                  }
                  tied.append(u);
               }
            }
         }
      }

      /*
      * Derivative of |k|^2 with respect to cell parameter n, averaged
      * over the tied minimum images in tied, or the derivative for
      * tied[0] if those of all images agree within a fraction 1.0E-9
      * (e.g., for the images k and -k of a wavevector that is half a
      * vector of the lattice of image shifts).
      */
      template <int D>
      double tiedImageDKSq(FSArray<IntVec<D>, 8> const & tied,
                           UnitCell<D> const & cell, int n)
      {
         const double epsilon = 1.0E-9;
         const int nTied = tied.size();
         const double value0 = cell.dksq(tied[0], n);
         double sum = value0;
         double maxDiff = 0.0;
         double scale = std::abs(value0);
         for (int t = 1; t < nTied; ++t) {
            double value = cell.dksq(tied[t], n);
            sum += value;
            maxDiff = std::max(maxDiff, std::abs(value - value0));
            scale = std::max(scale, std::abs(value));
         }
         if (maxDiff <= epsilon*scale) return value0;
         return sum/double(nTied);
      }

   }

   /*
   * Constructor.
   */
   template <int D>
   WaveList<D,CPT>::WaveList(bool isRealField)
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
   WaveList<D,CPT>::~WaveList() 
   {}

   /*
   * Allocate memory used by WaveList.
   */
   template <int D>
   void WaveList<D,CPT>::allocate(Mesh<D> const & m, UnitCell<D> const & c) 
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
         FFT<D,CPT>::computeKMesh(meshDimensions, kMeshDimensions_, kSize_);
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
               FFT<D,CPT>::hasImplicitInverse(kItr.position(), meshDimensions);
         }
      }

      clearUnitCellData();
      isAllocated_ = true;
   }

   /*
   * Clear all data that depends on unit cell parameters.
   */
   template <int D>
   void WaveList<D,CPT>::clearUnitCellData()
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
   void WaveList<D,CPT>::computeMinimumImages() 
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
   void WaveList<D,CPT>::computeKSq() 
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
   void WaveList<D,CPT>::computedKSq()
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
         RField<D,CPT>& dksq = dKSq_[i];
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
      
      // Minimum images tied with images whose |k|^2 has other derivatives
      // (cells with a variable angle only, e.g. rhombic at 120 degrees):
      // use the average over the tied images (see tiedImageDKSq)
      if (hasVariableAngle<D>(unitCell().lattice())) {
         FSArray<IntVec<D>, 14> neighbors;
         imageShiftNeighbors<D>(mesh().dimensions(), unitCell(), neighbors);
         FSArray<IntVec<D>, 8> tied;
         for (rank = 0; rank < kSize_; ++rank) {
            tiedMinimumImages<D>(minImages_[rank], neighbors, unitCell(),
                                 tied);
            if (tied.size() == 1) continue;
            for (i = 0; i < unitCell().nParameter(); ++i) {
               dKSq_[i][rank] = tiedImageDKSq<D>(tied, unitCell(), i);
               if (isRealField_) {
                  if (implicitInverse_[rank]) {
                     dKSq_[i][rank] *= 2.0;
                  }
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
   void WaveList<D,CPT>::sortWaves()
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
      int begin = 0;
      int end = 0;
      int ib, iw;
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
#endif
