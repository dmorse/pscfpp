#ifndef PSCF_MESH_ITERATOR_FORTRAN_TPP
#define PSCF_MESH_ITERATOR_FORTRAN_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MeshIteratorFortran.h"

namespace Pscf 
{

   template <int D>
   MeshIteratorFortran<D>::MeshIteratorFortran()
    : dimensions_(0),
      offsets_(0),
      position_(0),
      rank_(0),
      size_(0),
      atEnd_(false)
   {}

   template <int D>
   MeshIteratorFortran<D>::MeshIteratorFortran(IntVec<D> const& dimensions)
    : dimensions_(0),
      offsets_(0),
      position_(0),
      rank_(0),
      size_(0),
      atEnd_(false)
   {  setDimensions(dimensions); }

   template <int D>
   void MeshIteratorFortran<D>::setDimensions(IntVec<D> const& dimensions)
   {
      // Require that all dimensions are positive
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(dimensions[i] > 0);
      }
  
      // Set dimensions_ and size_ member variables
      dimensions_ = dimensions;
      size_ = 1;
      for (int i = 0; i < D; ++i) {
         size_ *= dimensions_[i];
      }

      // Compute offsets_ 
      offsets_[D - 1] = 1;
      for (int i = D - 1 ; i > 0; --i ) {
         offsets_[i - 1] = offsets_[i] * dimensions_[i];
      }

      // Initialize to first grid point
      begin();
   }

   template <int D>
   void MeshIteratorFortran<D>::begin()
   {
      UTIL_CHECK(size_ > 0);
      for (int i = 0; i < D; ++i) {
         position_[i] = 0;
      }
      rank_ = 0;
      atEnd_ = false;
   }

   template <int D>
   void MeshIteratorFortran<D>::operator++ () 
   {
      UTIL_CHECK(!atEnd_);

      // Cartesian index of position component that is incremented
      int k = 0;

      // Increment grid position vector
      bool next = true;
      while (next) {
         position_[k]++;
         if (position_[k] == dimensions_[k]) {
            position_[k] = 0;
            if (k == D - 1) {
               atEnd_ = true;
               next = false;
            } else {
               k++;
            }
         } else {
            next = false;
         }
      }

      // Compute array rank from updated position vector
      rank_ = 0;
      for (int dim = 0; dim < D; ++dim) {
         rank_ += offsets_[dim] * position_[dim];
      }

   }
   
}
#endif
