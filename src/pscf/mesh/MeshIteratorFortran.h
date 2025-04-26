#ifndef PSCF_MESH_ITERATOR_FORTRAN_H
#define PSCF_MESH_ITERATOR_FORTRAN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
   * and the integer rank of the current point in a mesh of specified
   * dimensions. The MeshIterator class instead iterates over points of
   * mesh in "C" order, in which the last index is most rapidly varying.
   *
   * Ranges of the position vector components and the scalar rank are 
   * defined using C/C++ conventions, exactly as in the Mesh and 
   * MeshIterator classes: Each component of the position vector is 
   * indexed from zero, and the scalar scalar rank of each mesh point is
   * defined using the C convention for the rank of a multi-dimensional 
   * array, in which the last index varies most rapidly. Because this
   * iterator visits points of mesh in Fortran order, in which the first
   * index varies most rapidly, the sequence of points visited by this 
   * iterator thus do not have sequential values for the rank.
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
      * Constructor, set mesh dimensions and initial iterator.
      *
      * Construction by this function is equivalent to construction by
      * the default constructor followed by a call to the setDimensions
      * function. On return, mesh dimensions are set and the iterator
      * is initialized to the first point in the mesh. 
      *
      * \param dimensions dimensions of the associated mesh
      */
      MeshIteratorFortran(IntVec<D> const& dimensions);

      /**
      * Set or reset the mesh dimensions, and initialize iterator.
      *
      * This fuction calls begin() internally. Upon return, the iterator
      * thus points to the first grid point, for which all position vector
      * components are equal to zero.
      *
      * \param dimensions dimensions of the associated mesh
      */
      void setDimensions(IntVec<D> const& dimensions);

      /**
      * Initialize the iterator to the first grid point.
      *
      * On return, the rank and all position components are all zero.
      */
      void begin();

      /**
      * Increment the iterator to the next grid point.
      */
      void operator ++ ();

      /**
      * Is this the end (i.e., past the last grid point)?
      */
      bool atEnd() const;

      /**
      * Return the scalar array rank of the associated grid point.
      */
      int rank() const;

      /**
      * Return a vector of coordinates of the associated grid point.
      *
      * Component i of the position varies from 0 to dimension[i] - 1,
      * for i = 0, ..., D - 1, where dimension[i] is the dimension of
      * the mesh in direction i. 
      */
      IntVec<D> const& position() const;

      /**
      * Return the mesh dimensions as an vector of integers.
      */
      IntVec<D> const& dimensions() const;

      /**
      * Return the mesh size (the number of grid points).
      */
      int size() const;

      /**
      * Return the vector of offsets.
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

