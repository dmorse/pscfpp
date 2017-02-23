#ifndef PSSP_MESH_ITERATOR_BASE_H
#define PSSP_MESH_ITERATOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>

namespace Pscf 
{

   using namespace Util;

   /**
   * Base class for mesh iterator class template.
   *
   * A mesh iterator iterates over the points of a mesh, keeping track
   * of both the IntVec<D> position and integer rank of the current
   * point as it goes. 
   *
   * \ingroup Pscf_Module
   */
   template <int D>
   class MeshIterator
   {

   public:

      /**
      * Default constructor
      */
      MeshIterator();

      /**
      * Constructor
      *
      * \param dimensions  IntVec<D> of grid dimensions
      */
      MeshIterator(const IntVec<D>& dimensions);

      // Compiler default destructor.

      // Compiler copy constructor.

      /**
      * Set the grid dimensions in all directions.
      *
      * \param dimensions  IntVec<D> of grid dimensions.
      */
      void setDimensions(const IntVec<D>& dimensions);

      /**
      * Set iterator to the first point in the mesh.
      */
      void begin();

      /**
      * Increment operator
      */
      void operator ++();

      /**
      * Get current position in the grid, as integer vector.
      */
      IntVec<D> position() const;

      /**
      * Get component i of the current position vector.
      *
      * \param i  index of Cartesian direction 0 <=i < D.
      */
      int position(int i) const;

      /**
      * Get the rank of current element.
      */
      int rank() const;

   private:

      /// Dimensions of grid
      IntVec<D> dimensions_;

      /// Current position in grid.
      IntVec<D> position_;

      /// Integer rank of current position.
      int rank_;

      /// Total number of grid points
      int size_;

      void increment(int i);

   };

   // Explicit specialization declarations

   template<> 
   void MeshIterator<1>::operator ++();

   template<> 
   void MeshIterator<2>::operator ++();

   template<> 
   void MeshIterator<3>::operator ++();

   // Inline member functions

   template <int D>
   inline IntVec<D> MeshIterator<D>::position() const
   {  return position_; }

   template <int D>
   inline int MeshIterator<D>::position(int i) const
   {
      assert(i >=0);
      assert(i < D);
      return dimensions_[i];
   }

   template <int D>
   inline int MeshIterator<D>::rank() const
   {  return rank_; }

   template <int D>
   inline void MeshIterator<D>::operator ++ ()
   {
      position_[D-1]++;
      if (position_[D-1] == dimensions_[D-1]) {
         position_[D-1] = 0;
         increment(D-2);
      }
      rank_++;
   }

   // Non-inline member functions

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
   {
      setDimensions(dimensions); 
   }

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

   template <int D>
   void MeshIterator<D>::increment(int i)
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
