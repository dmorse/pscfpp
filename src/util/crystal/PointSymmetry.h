#ifndef UTIL_POINT_SYMMETRY_H
#define UTIL_POINT_SYMMETRY_H

#include <util/space/Dimension.h>
#include <util/containers/FMatrix.h>

#include <iostream>

namespace Util {

   class IntVector;

   /**
   * A PointSymmetry represents a crystallographic point group symmetry.
   *
   * Crystallographic point group symmetries included 2, 3, and 4 fold rotations 
   * about axes, reflections through planes and inversion. Each such symmetry may 
   * be represented by a matrix of integers. 
   *
   * \ingroup Crystal_Module
   */
   class PointSymmetry 
   {

   public:

      /// Typedef for internal matrix
      typedef FMatrix<int, Dimension, Dimension> Matrix;

      /**
      * Default constructor.
      *
      * All elements of the rotation matrix are initialized to zero.
      */
      PointSymmetry();

      /**
      * Copy constructor.
      */
      PointSymmetry(const PointSymmetry& other);

      /**
      * Assignment operator.
      */
      PointSymmetry& operator = (const PointSymmetry& other);

      /**
      * Return the inverse of this symmetry element.
      */
      PointSymmetry inverse() const;

      /**
      * Return an element of the matrix by reference.
      */
      int& R(int i, int j);
 
      /**
      * Return an element of the matrix by value
      */
      int  R(int i, int j) const;

      // Static method

      /**
      * Return the identity element.
      */
      static const PointSymmetry& identity();
 
   private:

      /**
      * Rotation matrix. 
      *
      * For a point group symmetry, elements of this matrix are all 0, 1 or -1.
      * There should be only one nonzero element per row or column, so that the
      * determinant is 1 or -1. 
      */
      Matrix R_;

      /// Identity element (static member stored for reference)
      static PointSymmetry identity_;

      /// Has the static identity_ been constructed?
      static bool hasIdentity_;

      /// Construct static identity_ object.
      static void makeIdentity();

   // friends:

      friend bool operator==(const PointSymmetry& A, const PointSymmetry& B);

      friend bool operator!=(const PointSymmetry& A, const PointSymmetry& B);

      friend PointSymmetry 
      operator * (const PointSymmetry& A, const PointSymmetry& B);

      friend IntVector operator * (const PointSymmetry& S, const IntVector& V);

      friend IntVector operator * (const IntVector& V, const PointSymmetry& S);

      friend std::ostream& operator << (std::ostream& out, const PointSymmetry& A);

   };

   // Friend function declarations

   /**
   * Are two PointSymmetry objects equivalent?
   *
   * \param A First  object
   * \param B Second object
   * \return True if A == B, false otherwise
   */
   bool operator == (const PointSymmetry& A, const PointSymmetry& B);

   /**
   * Are two PointSymmetry objects not equivalent?
   *
   * \param A First  object
   * \param B Second object
   * \return True if A != B, false otherwise
   */
   bool operator != (const PointSymmetry& A, const PointSymmetry& B);

   /**
   * Return the product A*B of two symmetry objects.
   *
   * \param A First  object
   * \param B Second object
   * \return product A*B
   */
   PointSymmetry operator * (const PointSymmetry& A, const PointSymmetry& B);

   /**
   * Return the IntVector product S*V of a rotation matrix and an IntVector.
   *
   * \param S PointSymmetry
   * \param V Integer Vector
   * \return product S*V
   */
   IntVector operator * (const PointSymmetry& S, const IntVector& V);

   /**
   * Return the IntVector product V*S of an IntVector and a rotation matrix.
   *
   * \param S PointSymmetry
   * \param V Integer Vector
   * \return product S*V
   */
   IntVector operator * (const IntVector& V, const PointSymmetry& S);

   /**
   * Output stream inserter for a PointSymmetry.
   *
   * \param out output stream
   * \param A   PointSymmetry object to be output
   * \return modified output stream
   */ 
   std::ostream& operator << (std::ostream& out, const PointSymmetry& A);

   // Inline method definitions

   inline bool operator != (const PointSymmetry& A, const PointSymmetry& B)
   {  return !(A == B); }

   /*
   * Return an element of the matrix by reference.
   */
   inline int& PointSymmetry::R(int i, int j)
   {  return R_(i, j); }
 
   /*
   * Return an element of the matrix by value
   */
   inline int  PointSymmetry::R(int i, int j) const
   {  return R_(i, j); }

   /*
   * Return identity element.
   */
   inline const PointSymmetry& PointSymmetry::identity()
   {
      if (!hasIdentity_) makeIdentity();
      return identity_;
   }
 
}
#endif
