#ifndef PSCF_SPACE_SYMMETRY_H
#define PSCF_SPACE_SYMMETRY_H

#include <pscf/math/IntVec.h>
#include <util/math/Rational.h>
#include <util/containers/FMatrix.h>
#include <util/containers/FArray.h>

#include <iostream>

namespace Pscf {

   using namespace Util;

   template <int D> class SpaceSymmetry;

   // Non-member function declarations

   /**
   * Are two SpaceSymmetry objects equivalent?
   *
   * \param A first symmetry
   * \param B second symmetry
   * \return True if A == B, false otherwise
   */
   template <int D>
   bool operator == (const SpaceSymmetry<D>& A, const SpaceSymmetry<D>& B);

   /**
   * Are two SpaceSymmetry objects not equivalent?
   *
   * \param A first symmetry
   * \param B second symmetry
   * \return True if A != B, false otherwise
   */
   template <int D>
   bool operator != (const SpaceSymmetry<D>& A, const SpaceSymmetry<D>& B);

   /**
   * Return the product A*B of two symmetry objects.
   *
   * \param A first symmetry
   * \param B second symmetry
   * \return product A*B
   */
   template <int D>
   SpaceSymmetry<D> 
   operator * (const SpaceSymmetry<D>& A, const SpaceSymmetry<D>& B);

   /**
   * Return the IntVec<D> product S*V of a rotation matrix and an IntVec<D>.
   * 
   * The product is defined to be the matrix product of the rotation matrix
   * and the integer vector S.R * V.
   *
   * \param S symmetry operation
   * \param V integer vector
   * \return product S*V
   */
   template <int D>
   IntVec<D> operator * (const SpaceSymmetry<D>& S, const IntVec<D>& V);

   /**
   * Return the IntVec<D> product V*S of an IntVec<D> and a rotation matrix.
   *
   * The product is defined to be the matrix product of the integer vector
   * and the space group rotation matrix S.R * V.
   *
   * \param V integer vector
   * \param S symmetry operation
   * \return product V*S
   */
   template <int D>
   IntVec<D> operator * (const IntVec<D>& V, const SpaceSymmetry<D>& S);

   /**
   * Output stream inserter for a SpaceSymmetry<D>
   *
   * \param out output stream
   * \param A  SpaceSymmetry<D> object to be output
   * \return  modified output stream
   */ 
   template <int D>
   std::ostream& operator << (std::ostream& out, const SpaceSymmetry<D>& A);

   /**
   * Input stream extractor for a SpaceSymmetry<D>
   *
   * \param in  input stream
   * \param A  SpaceSymmetry<D> object to be input
   * \return  modified input stream
   */ 
   template <int D>
   std::istream& operator >> (std::istream& in, SpaceSymmetry<D>& A);

   /**
   * A SpaceSymmetry represents a crystallographic space group symmetry.
   *
   * Crystallographic space group symmetry operation combines a point group 
   * operation (e.g., 2, 3, and 4 fold rotations about axes, reflections, or
   * inversion) with possible translations by a fraction of a unit cell.
   *
   * Both the rotation matrix R and the translation t are represented using 
   * a basis of Bravais lattice basis vectors.  Because Bravais basis vectors 
   * must map onto other lattice vectors, this implies that elements of all 
   * elements of the rotation matrix must be integers.  To guarantee that the 
   * inverse of the rotation matrix is also a matrix of integers, we require 
   * that the determinant of the rotation matrix must be +1 or -1. The
   * translation vector is represented by a vector of D rational numbers
   * (i.e., fractions) of the form n/m with m = 2, 3, or 4 and n < m. 
   *
   * The basis used to describe a crytallographic group may be either a
   * primitive or non-primitive unit cell. Thus, for example, the space
   * group of a bcc crystal may be expressed either using a basis of 3
   * three orthogonal simple cubic unit vectors, with a translation
   * t = (1/2, 1/2, 1/2), or as a point group using a set of three 
   * non-orthogonal basis vectors for the primitive unit cell. 
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   class SpaceSymmetry 
   {

   public:

      /// Typedef for matrix used to represent point group operation.
      typedef FMatrix<int, D, D> Rotation;

      /// Typedef for vector used to represent fractional translation.
      typedef FArray<Rational, D> Translation;

      /**
      * Default constructor.
      *
      * All elements of the rotation matrix are initialized to zero.
      */
      SpaceSymmetry();

      /**
      * Copy constructor.
      */
      SpaceSymmetry(const SpaceSymmetry<D>& other);

      /**
      * Assignment operator.
      */
      SpaceSymmetry<D>& operator = (const SpaceSymmetry<D>& other);

      /**
      * Shift components of translation to [0,1).
      */
      void normalize();

      /**
      * Compute and return the inverse of this symmetry element.
      */
      SpaceSymmetry<D> inverse() const;

      /**
      * Compute and return the inverse of the rotation matrix.
      */
      SpaceSymmetry<D>::Rotation inverseRotation() const;

      /**
      * Compute and return the determinant of the rotation matrix.
      */
      int determinant() const;

      /**
      * Return an element of the matrix by reference.
      *
      * \param i  1st (row) index
      * \param j  2nd (column) index
      */
      int& R(int i, int j);
 
      /**
      * Return an element of the matrix by value.
      *
      * \param i  1st (row) index
      * \param j  2nd (column) index
      */
      int R(int i, int j) const;

      /**
      * Return a component of the translation by reference.
      * 
      * \param i component index 
      */
      Rational& t(int i);
 
      /**
      * Return an element of the translation by value.
      *
      * \param i component index 
      */
      Rational t(int i) const;

      // Static method

      /**
      * Return the identity element.
      */
      static const SpaceSymmetry<D>& identity();

   private:

      /**
      * Matrix representation of point group operation.
      *
      * Because it is expressed in a Bravais basis, and Bravais lattice 
      * vectors must map onto other Bravais lattice vectors, elements of 
      * this matrix are integers. 
      */
      Rotation R_;

      /**
      * Translation vector
      */
      Translation t_;

      // Static member variables

      /// Identity element (static member stored for reference)
      static SpaceSymmetry<D> identity_;

      /// Has the static identity_ been constructed?
      static bool hasIdentity_;

      /// Construct static identity_ object.
      static void makeIdentity();

   // friends:

      friend 
      bool operator == <> (const SpaceSymmetry<D>& A, 
                           const SpaceSymmetry<D>& B);

      friend 
      bool operator != <> (const SpaceSymmetry<D>& A, 
                           const SpaceSymmetry<D>& B);

      friend 
      SpaceSymmetry<D>
      operator * <> (const SpaceSymmetry<D>& A, 
                     const SpaceSymmetry<D>& B);

      friend 
      IntVec<D> operator * <> (const IntVec<D>& V, 
                               const SpaceSymmetry<D>& S);

      friend 
      IntVec<D> operator * <>(const SpaceSymmetry<D>& S, const IntVec<D>& V);

      friend 
      std::ostream& operator << <> (std::ostream& out, 
                                    const SpaceSymmetry<D>& A);

      friend 
      std::istream& operator >> <> (std::istream& in, 
                                    SpaceSymmetry<D>& A);

   };

   // Static member variable declaration templates

   template <int D> 
   SpaceSymmetry<D> SpaceSymmetry<D>::identity_;

   template <int D> 
   bool SpaceSymmetry<D>::hasIdentity_;

   // Explicit specialization of members

   template <> 
   SpaceSymmetry<1>::Rotation SpaceSymmetry<1>::inverseRotation() const;

   template <> 
   SpaceSymmetry<2>::Rotation SpaceSymmetry<2>::inverseRotation() const;

   template <> 
   SpaceSymmetry<3>::Rotation SpaceSymmetry<3>::inverseRotation() const;

   template <> 
   int SpaceSymmetry<1>::determinant() const;

   template <> 
   int SpaceSymmetry<2>::determinant() const;

   template <> 
   int SpaceSymmetry<3>::determinant() const;

   // Inline method definitions

   template <int D>
   inline 
   bool operator != (const SpaceSymmetry<D>& A, const SpaceSymmetry<D>& B)
   {  return !(A == B); }

   /*
   * Return an element of the matrix by reference.
   */
   template <int D>
   inline 
   int& SpaceSymmetry<D>::R(int i, int j)
   {  return R_(i, j); }
 
   /*
   * Return an element of the matrix by value
   */
   template <int D>
   inline 
   int SpaceSymmetry<D>::R(int i, int j) const
   {  return R_(i, j); }

   /*
   * Return an element of the translation vector by reference.
   */
   template <int D>
   inline 
   Rational& SpaceSymmetry<D>::t(int i)
   {  return t_[i]; }
 
   /*
   * Return an element of the translation vector by value
   */
   template <int D>
   inline 
   Rational SpaceSymmetry<D>::t(int i) const
   {  return t_[i]; }

   /*
   * Return the identity symmetry operation.
   */
   template <int D>
   inline 
   const SpaceSymmetry<D>& SpaceSymmetry<D>::identity()
   {
      if (!hasIdentity_) makeIdentity();
      return identity_;
   }
 
}
#include "SpaceSymmetry.tpp"
#endif
