#ifndef PSCF_MESH_ITERATOR_FORTRAN_H
#define PSCF_MESH_ITERATOR_FORTRAN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>   // member

namespace Pscf 
{

   using namespace Util;

   /**
   * Iterator over points in a mesh in "Fortran" order.
   *
   * This mesh iterator iterates over the points of a mesh in "Fortran"
   * order, in which the first array index varies most rapidly. During 
   * iteration, the class keeps track of both the IntVec<D> position 
   * and the integer rank of that point in the mesh. 
   *
   * Ranges of the position components and the rank are defined using
   * C conventions, exactly as in the Mesh and MeshIterator classes. 
   * Each component of the position vector is indexed from zero, and the 
   * scalar rank of each mesh point is defined using the C convention 
   * for the rank of a multi-dimensional array, in which the last index 
   * varies most rapidly.  The sequence of points visited by this iterator
   * thus do not have sequential values for the rank.
   *
   * Typical usage:
   * \code
   *    MeshIteratorFortran<D> iter(meshDimensions);
   *    for (iter.begin(); !iter.atEnd(); ++iter()) {
   *       (Do something) 
   *       cout << iter.position() << iter.rank();
   *    }
   * \endcode
   *
   * \ingroup Pscf_Mesh_Module
   */
   template <int D>
   class MeshIteratorFortran 
   {

   public:

      /**
      * Default constructor.
      */
      MeshIteratorFortran();

      /**
      * Constructor, set mesh dimensions
      *
      * \param dimensions dimensions of the associated mesh
      */
      MeshIteratorFortran(IntVec<D> const& dimensions);

      /**
      * Set or reset the mesh dimensions, and initialize iterator.
      *
      * This fuction calls begin() internally. Upon return the
      * iterator thus points to the first grid point, with all
      * position coordinates equal to zero.
      *
      * \param dimensions dimensions of the associated mesh
      */
      void setDimensions(IntVec<D> const& dimensions);

      /**
      * Initialize the iterator to the first grid point.
      *
      * On return, rank and position components are all zero.
      */
      void begin();

      /**
      * Increment the position and rank to the next grid point.
      */
      void operator ++ ();

      /**
      * Is this the end (i.e., past the last grid point)?
      */
      bool atEnd() const;

      /**
      * Return the array rank of the associated grid point.
      */
      int rank() const;

      /**
      * Return integer coordinates of the associated grid point.
      */
      IntVec<D> const& position() const;

      /**
      * Return mesh dimensions.
      */
      IntVec<D> const& dimensions() const;

      /**
      * Return the mesh size (the number of grid points).
      */
      int size() const;

      /**
      * Return vector of offsets.
      */
      IntVec<D> const& offsets() const;

   private:

      // Dimensions of associated Mesh<D> object
      IntVec<D> dimensions_;
 
      // Offsets associated with position components
      IntVec<D> offsets_;

      // Grid position - vector of integer components 
      IntVec<D> position_;

      // Scalar rank of current node
      int  rank_;

      // Mesh size
      int  size_;

      // Has the iterator past the last grid point?
      bool atEnd_;
      
   };
   
   template <int D>
   inline bool MeshIteratorFortran<D>::atEnd() const
   {  return atEnd_; }

   template <int D>
   inline int MeshIteratorFortran<D>::rank() const
   {  return rank_; }
   
   template <int D>
   inline IntVec<D> const& MeshIteratorFortran<D>::position() const
   {  return position_; }

   template <int D>
   inline IntVec<D> const& MeshIteratorFortran<D>::dimensions() const
   {  return dimensions_; }

   template <int D>
   inline int MeshIteratorFortran<D>::size() const
   {  return size_; }
   
   template <int D>
   inline IntVec<D> const& MeshIteratorFortran<D>::offsets() const
   {  return offsets_; }

   #ifndef PSCF_MESH_ITERATOR_FORTRAN_TPP
   // Suppress implicit instantiation
   extern template class MeshIteratorFortran<1>;
   extern template class MeshIteratorFortran<2>;
   extern template class MeshIteratorFortran<3>;
   #endif

}
#endif

