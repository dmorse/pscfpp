#ifndef PSCF_MESH_ITERATOR_H
#define PSCF_MESH_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>           // member

namespace Pscf 
{

   using namespace Util;

   /**
   * Iterator over points in a Mesh<D>.
   *
   * A mesh iterator iterates over the points of a mesh, keeping track
   * of both the IntVec<D> position and integer rank of the current
   * point as it goes. 
   *
   * \ingroup Pscf_Mesh_Module
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

      // Compiler copy constructor.

      // Compiler default destructor.

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
      * Increment iterator to next mesh point.
      */
      void operator ++();

      /**
      * Is this the end (i.e., one past the last point)?
      */
      bool atEnd() const;

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

      /// Recursive function for multi-dimensional increment.
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

   // Return the entire IntVec<D>
   template <int D>
   inline IntVec<D> MeshIterator<D>::position() const
   {  return position_; }

   // Return one component of the position IntVec<D>
   template <int D>
   inline int MeshIterator<D>::position(int i) const
   {
      assert(i >=0);
      assert(i < D);
      return position_[i];
   }

   // Return the rank
   template <int D>
   inline int MeshIterator<D>::rank() const
   {  return rank_; }

   // Is this the end (i.e., one past the last point)?
   template <int D>
   inline bool MeshIterator<D>::atEnd() const
   { return (bool)(rank_ == size_); }

   // Inline explicit specializations 

   template <>
   inline void MeshIterator<1>::operator ++ ()
   {
      position_[0]++;
      if (position_[0] == dimensions_[0]) {
         position_[0] = 0;
      }
      rank_++;
   }

   template <>
   inline void MeshIterator<2>::operator ++ ()
   {
      position_[1]++;
      if (position_[1] == dimensions_[1]) {
         position_[1] = 0;
         increment(0);
      }
      rank_++;
   }

   template <>
   inline void MeshIterator<3>::operator ++ ()
   {
      position_[2]++;
      if (position_[2] == dimensions_[2]) {
         position_[2] = 0;
         increment(1);
      }
      rank_++;
   }

   #ifndef PSCF_MESH_ITERATOR_TPP
   // Suppress implicit instantiation
   extern template class MeshIterator<1>;
   extern template class MeshIterator<2>;
   extern template class MeshIterator<3>;
   #endif

}
#endif
