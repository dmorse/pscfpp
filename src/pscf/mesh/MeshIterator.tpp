#ifndef PSCF_MESH_ITERATOR_TPP
#define PSCF_MESH_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MeshIterator.h"

namespace Pscf 
{

   /*
   * Default constructor
   */ 
   template <int D>
   MeshIterator<D>::MeshIterator()
    : dimensions_(0),
      position_(0),
      rank_(0),
      size_(0)
   {}

   /*
   * Constructor
   */ 
   template <int D>
   MeshIterator<D>::MeshIterator(const IntVec<D>& dimensions)
    : dimensions_(0),
      position_(0),
      rank_(0),
      size_(0)
   {  setDimensions(dimensions); }

   /*
   * Set the mesh dimensions.
   */
   template <int D>
   void MeshIterator<D>::setDimensions(const IntVec<D>& dimensions)
   { 
      for (int i = 0; i < D; ++i) {
         if (dimensions[i] <= 0) {
            UTIL_THROW("Mesh dimensions must be positive");
         }
      }
 
      dimensions_ = dimensions;
      size_ = 1;
      for (int i = 0; i < D; ++i) {
         size_ *= dimensions_[i];
      }
   }

   /*
   * Reset iterator to point to first element.
   */
   template <int D>
   void MeshIterator<D>::begin()
   {
      rank_ = 0;
      for (int i = 0; i < D; ++i) {
         position_[i] = 0;
      }
   }

   /*
   * Increment this iterator to the next mesh point (default templ.).
   *
   * Note: Explicit instantiations are defined for D = 1, 2 and 3.
   * This default implementation should thus not normally be used.
   */
   template <int D>
   inline void MeshIterator<D>::operator ++()
   {
      position_[D-1]++;
      if (position_[D-1] == dimensions_[D-1]) {
         position_[D-1] = 0;
         if (D > 1) {
            increment(D-2);
         }
      }
      rank_++;
   }

   /*
   * Recursive increment function (private).
   */ 
   template <int D>
   inline void MeshIterator<D>::increment(int i)
   {
      position_[i]++;
      if (position_[i] == dimensions_[i]) {
         position_[i] = 0;
         if (i > 0) {
            increment(i-1);
         }
      }
   }

}
#endif
