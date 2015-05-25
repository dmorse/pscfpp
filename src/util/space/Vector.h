#ifndef UTIL_VECTOR_H
#define UTIL_VECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/space/Dimension.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
#endif

#include <iostream>
#include <cmath>

namespace Util
{

   /**
   * A Vector is a Cartesian vector.
   *
   * The Cartesian elements of a Vector can be accessed using array notation:
   * The elements of a three dimensional Vector v are v[0], v[1], and v[2].
   * The subscript operator [] returns elements as references, which can be
   * used on either the left or right side of an assignment operator.
   *
   * The arithmetic assignment operators +=, -=, *=, and /= are overloaded.
   * The operators += and -= represent increment or decrement by a vector,
   * while *= and /= represent multiplication or division by a scalar.
   *
   * All other unary and binary mathematical operations are implemented as
   * methods. Operations that yield a scalar result, such as a dot product,
   * are implemented as methods that return the resulting value. Operations
   * that yield a Vector, such as vector addition, are implemented by methods
   * that assign the result to the invoking vector, and return a reference
   * to the invoking vector. For example,
   * \code
   *    Vector a, b, c;
   *    double s;
   *
   *    a[0] = 0.0
   *    a[1] = 1.0
   *    a[2] = 2.0
   *
   *    b[0] =  0.5
   *    b[1] = -0.5
   *    b[2] = -1.5
   *
   *    // Set s = a.b
   *    s = a.dot(b)
   *
   *    // Set c = a + b
   *    c.add(a, b)
   *
   *    // Set a = a + b
   *    a += b
   *
   *    // Set b = b*2
   *    b *= 2
   *
   * \endcode
   * This syntax for Vector valued operations avoids dynamic allocation of
   * temporary Vector objects, by requiring that the invoking function
   * provide an object to hold the result.
   *
   * For efficiency, all methods in this class are declared inline.
   *
   * \ingroup Space_Module
   */
   class Vector
   {

   public:

      /// \name Constructors
      //@{

      /**
      * Default constructor
      */
      Vector();

      /**
      * Copy constructor
      *
      * \param v Vector to be copied
      */
      Vector(const Vector& v);

      /**
      * Constructor, initialize all elements to a scalar value.
      *
      * \param scalar initial value for all elements.
      */
      explicit Vector(double scalar);

      /**
      * Construct Vector from C double[3] array.
      *
      * \param v array of 3 coordinates
      */
      explicit Vector(const double* v);

      /**
      * Construct Vector from its coordinates.
      *
      * \param x x-axis coordinate, v[0]
      * \param y y-axis coordinate, v[1]
      * \param z z-axis coordinate, v[2]
      */
      Vector(double x, double y, double z=0.0);

      //@}

      /**
      * Set all elements of a 3D vector to zero.
      */
      Vector& zero();

      /**
      * Serialize to/from an archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /// \name Assignment
      //@{

      /**
      * Copy assignment.
      *
      * \param v Vector to assign.
      */
      Vector& operator = (const Vector& v);

      /**
      * Assignment from C double[3] array.
      *
      * \param v  C-array of components
      */
      Vector& operator = (const double* v);

      //@}
      /// \name Arithmetic Assignment
      //@{

      /**
      * Add vector dv to this vector.
      *
      * Upon return, *this = this + dv.
      *
      * \param dv  vector increment (input)
      */
      void operator += (const Vector& dv);

      /**
      * Subtract vector dv from this vector.
      *
      * Upon return, *this = this + dv.
      *
      * \param dv  vector increment (input)
      */
      void operator -= (const Vector& dv);

      /**
      * Multiply this vector by scalar s.
      *
      * Upon return, *this = (*this)*s.
      *
      * \param s  scalar multiplier
      */
      void operator *= (double s);

      /**
      * Divide this vector by scalar s.
      *
      * Upon return, *this = (*this)/s.
      *
      * \param s  scalar divisor (input)
      */
      void operator /= (double s);

      //@}
      /// \name Array Subscript
      //@{

      /**
      * Return one Cartesian element by value.
      *
      * \param i  element index
      * \return element i of the vector
      */
      const double& operator [] (int i) const;

      /**
      * Return one element of the vector by references.
      *
      * \param i  element index
      * \return element i of this vector
      */
      double& operator [] (int i);

      //@}
      /// \name Scalar-valued functions
      //@{

      /**
      * Return square magnitude of this vector.
      *
      * \return square magnitude of this vector
      */
      double square() const;

      /**
      * Return absolute magnitude of this vector.
      *
      * \return absolute magnitude (norm) of this vector.
      */
      double abs() const;

      /**
      * Return dot product of this vector and vector v.
      *
      * \param  v input vector
      * \return dot product of this vector and vector v
      */
      double dot(const Vector& v) const;

      /**
      * Return projection of this vector along vector p.
      *
      * \param  p vector parallel to direction along which to project
      * \return scalar projection this->dot(p)/p.abs()
      */
      double projection(const Vector& p) const;

      //@}
      /// \name Vector valued functions (result assigned to invoking object)
      //@{

      /**
      * Add vectors v1 and v2.
      *
      * Upon return, *this = v1 + v2.
      *
      * \param v1  vector (input)
      * \param v2  vector (input)
      * \return modified invoking vector
      */
      Vector& add(const Vector& v1, const Vector& v2);

      /**
      * Subtract vector v2 from v1.
      *
      * Upon return, *this = v1 - v2.
      *
      * \param v1  vector (input)
      * \param v2  vector (input)
      * \return modified invoking vector
      */
      Vector& subtract(const Vector& v1, const Vector& v2);

      /**
      * Multiply a vector v by a scalar s.
      *
      * Upon return, *this = v*s.
      *
      * \param v  vector input
      * \param s  scalar input
      * \return modified invoking vector
      */
      Vector& multiply(const Vector& v, double s);

      /**
      * Divide vector v by scalar s.
      *
      * Upon return, *this = v/s;
      *
      * \param v  vector input
      * \param s  scalar input
      * \return modified invoking vector
      */
      Vector& divide(const Vector& v, double s);

      /**
      * Calculate cross product of vectors v1 and v2.
      *
      * Upon return, *this = v1 x v2.
      *
      * \param v1  input vector
      * \param v2  input vector
      * \return modified invoking vector
      */
      Vector& cross(const Vector& v1, const Vector& v2);

      /**
      * Calculate unit vector parallel to input vector v.
      *
      * Upon return *this = unit vector.
      *
      * \param v input vector
      * \return modified invoking Vector
      */
      Vector& versor(const Vector& v);

      /**
      * Calculate component of vector v parallel to vector p.
      *
      * Upon return, the invoking vector is equal to the vector projection of
      * vector v along a direction parallel to vector p.
      *
      * The vector projection of v along p is parallel to p and has an absolute
      * magnitude equal to the scalar projection of v along p.
      *
      * \param v vector to project
      * \param p vector along which to project
      * \return modified invoking Vector
      */
      Vector& parallel(const Vector& v, const Vector& p);

      /**
      * Calculate component of vector v transverse to vector p.
      *
      * Upon return, the invoking vector is equal to the vector
      * projection of vector v perpendicular to vector p.
      *
      * \param v  input vector
      * \param p  vector perpendicular to which to project.
      * \return modified invoking Vector
      */
      Vector& transverse(const Vector& v, const Vector& p);

      /**
      * Computes the index corresponding to minimum element in a vector.
      *
      * \param v  input vector
      * \return index of the minimum element.
      */
      int minId(const Vector& v);

      /**
      * Computes the index corresponding to maximum element in a vector.
      *
      * \param v  input vector
      * \return index of the maximum element.
      */
      int maxId(const Vector& v);

      //@}
      /// \name Static Members
      //@{

      /**
      * Initialize Zero Vector.
      */
      static void initStatic();

      /**
      * Zero Vector = {0.0, 0.0, 0.0}
      */
      static const Vector Zero;

      #ifdef UTIL_MPI
      /**
      * Commit MPI datatype MpiTraits<Vector>::type.
      */
      static void commitMpiType();
      #endif

      //@}

   private:

      /// Width of field per Cartesian coordinate in stream IO
      static const int Width = 25;

      /// Precision in stream IO of Vector coordinates
      static const int Precision = 17;

      /// Elements of the vector.
      double elem_[Dimension];

   //friends:

      friend bool operator == (const Vector& v1, const Vector& v2);

      friend bool operator == (const Vector& v1, const double* v2);

      friend 
      std::istream& operator >> (std::istream& in, Vector &vector);

      friend 
      std::ostream& operator << (std::ostream& out, const Vector &vector);

   };

   /// Equality for Vectors.
   bool operator == (const Vector& v1, const Vector& v2);

   /// Equality of Vector and C array.
   bool operator == (const Vector& v1, const double* v2);

   /// Equality of C array and Vector.
   bool operator == (const double* v1, const Vector& v2);

   /// Inequality of two Vectors.
   bool operator != (const Vector& v1, const Vector& v2);

   /// Inequality of Vector and C array.
   bool operator != (const Vector& v1, const double* v2);

   /// Inequality of C array and Vector.
   bool operator != (const double* v1, const Vector& v2);

   /**
   * istream extractor for a Vector.
   *
   * Input elements of a vector from stream, without line breaks.
   *
   * \param in      input stream
   * \param vector  Vector to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Vector &vector);

   /**
   * ostream inserter for a Vector.
   *
   * Output elements of a vector to stream, without line breaks.
   * \param  out     output stream
   * \param  vector  Vector to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Vector &vector);

   #ifdef UTIL_MPI
   /**
   * Explicit specialization MpiTraits<Vector>.
   */
   template <>
   class MpiTraits<Vector>
   {
   public:
      static MPI::Datatype type;   ///< MPI Datatype
      static bool hasType;         ///< Is the MPI type initialized?
   };
   #endif

   // Inline methods

   /*
   * Default constructor
   */
   inline
   Vector::Vector()
   {}

   /*
   * Copy constructor
   */
   inline
   Vector::Vector(const Vector& v)
   {
      elem_[0] = v.elem_[0];
      elem_[1] = v.elem_[1];
      elem_[2] = v.elem_[2];
   }

   /*
   * Constructor, initialize all elements to a scalar value.
   */
   inline
   Vector::Vector(double scalar)
   {
      elem_[0] = scalar;
      elem_[1] = scalar;
      elem_[2] = scalar;
   }

   /*
   * Construct Vector from C double[3] array.
   */
   inline
   Vector::Vector(const double* v)
   {
      elem_[0] = v[0];
      elem_[1] = v[1];
      elem_[2] = v[2];
   }

   /*
   * Construct Vector from its coordinates.
   */
   inline
   Vector::Vector(double x, double y, double z)
   {
      elem_[0] = x;
      elem_[1] = y;
      elem_[2] = z;
   }

   /*
   * Set all elements of a 3D vector to zero.
   */
   inline
   Vector& Vector::zero()
   {
      elem_[0] = 0.0;
      elem_[1] = 0.0;
      elem_[2] = 0.0;
      return *this;
   }

   /*
   * Assignment.
   */
   inline
   Vector& Vector::operator = (const Vector& v)
   {
      elem_[0] = v.elem_[0];
      elem_[1] = v.elem_[1];
      elem_[2] = v.elem_[2];
      return *this;
   }

   /*
   * Assignment from C double[] array.
   */
   inline
   Vector& Vector::operator = (const double* v)
   {
      elem_[0] = v[0];
      elem_[1] = v[1];
      elem_[2] = v[2];
      return *this;
   }

   /*
   * Add vector dv to this vector.
   */
   inline
   void Vector::operator += (const Vector& dv)
   {
      elem_[0] += dv.elem_[0];
      elem_[1] += dv.elem_[1];
      elem_[2] += dv.elem_[2];
   }

   /*
   * Subtract vector dv from this vector.
   */
   inline
   void Vector::operator -= (const Vector& dv)
   {
      elem_[0] -= dv.elem_[0];
      elem_[1] -= dv.elem_[1];
      elem_[2] -= dv.elem_[2];
   }

   /*
   * Multiply this vector by scalar s.
   */
   inline
   void Vector::operator *= (double s)
   {
      elem_[0] *= s;
      elem_[1] *= s;
      elem_[2] *= s;
   }

   /*
   * Divide this vector by scalar s.
   */
   inline
   void Vector::operator /= (double s)
   {
      elem_[0] /= s;
      elem_[1] /= s;
      elem_[2] /= s;
   }

   /*
   * Return one Cartesian element by value.
   */
   inline
   const double& Vector::operator [] (int i) const
   { 
      assert(i >=0); 
      assert(i < Dimension); 
      return elem_[i]; 
   }

   /*
   * Return a reference to one element of the vector.
   */
   inline
   double& Vector::operator [] (int i)
   {
      assert(i >=0); 
      assert(i < Dimension); 
      return elem_[i]; 
   }

   /*
   * Return square magnitude of this vector.
   */
   inline
   double Vector::square() const
   {
      return (elem_[0]*elem_[0] + elem_[1]*elem_[1] + elem_[2]*elem_[2]);
   }

   /*
   * Return absolute magnitude of this vector.
   */
   inline
   double Vector::abs() const
   {  return sqrt(square()); }

   /*
   * Return dot product of this vector and vector v.
   */
   inline
   double Vector::dot(const Vector& v) const
   {
      return elem_[0]*v.elem_[0] + elem_[1]*v.elem_[1] + elem_[2]*v.elem_[2];
   }

   /*
   * Return projection of this vector along vector p.
   */
   inline
   double Vector::projection(const Vector& p) const
   {
      double abs_p = p.abs();
      if (abs_p > 1.0E-8) {
         return dot(p)/abs_p;
      } else {
         return 0.0;
      }
   }

   /*
   * Add vectors v1 and v2.
   *
   * Upon return, *this = v1 + v2.
   */
   inline
   Vector& Vector::add(const Vector& v1, const Vector& v2)
   {
      elem_[0] = v1.elem_[0] + v2.elem_[0];
      elem_[1] = v1.elem_[1] + v2.elem_[1];
      elem_[2] = v1.elem_[2] + v2.elem_[2];
      return *this;
   }

   /*
   * Subtract vector v2 from v1.
   *
   * Upon return, *this = v1 - v2.
   * \return modified invoking vector
   */
   inline
   Vector& Vector::subtract(const Vector& v1, const Vector& v2)
   {
      elem_[0] = v1.elem_[0] - v2.elem_[0];
      elem_[1] = v1.elem_[1] - v2.elem_[1];
      elem_[2] = v1.elem_[2] - v2.elem_[2];
      return *this;
   }

   /*
   * Multiply a vector v by a scalar s.
   *
   * Upon return, *this = v*s.
   */
   inline
   Vector& Vector::multiply(const Vector& v, double s)
   {
      elem_[0] = v.elem_[0]*s;
      elem_[1] = v.elem_[1]*s;
      elem_[2] = v.elem_[2]*s;
      return *this;
   }

   /*
   * Divide vector v by scalar s.
   *
   * Upon return, *this = v/s;
   */
   inline
   Vector& Vector::divide(const Vector& v, double s)
   {
      elem_[0] = v.elem_[0]/s;
      elem_[1] = v.elem_[1]/s;
      elem_[2] = v.elem_[2]/s;
      return *this;
   }

   /*
   * Calculate cross product of vectors v1 and v2.
   *
   * Upon return, *this = v1 x v2.
   */
   inline
   Vector& Vector::cross(const Vector& v1, const Vector& v2)
   {
      elem_[0] = v1.elem_[1]*v2.elem_[2] - v1.elem_[2]*v2.elem_[1];
      elem_[1] = v1.elem_[2]*v2.elem_[0] - v1.elem_[0]*v2.elem_[2];
      elem_[2] = v1.elem_[0]*v2.elem_[1] - v1.elem_[1]*v2.elem_[0];
      return *this;
   }

   /*
   * Calculate unit vector parallel to v: this = v/|v|
   */
   inline
   Vector& Vector::versor(const Vector& v)
   {
      double vabs = v.abs();
      if (vabs > 1.0E-8) {
         vabs = 1.0/vabs;
         elem_[0] = v.elem_[0]*vabs;
         elem_[1] = v.elem_[1]*vabs;
         elem_[2] = v.elem_[2]*vabs;
      } else {
         elem_[0] = 1.0;
         elem_[1] = 0.0;
         elem_[2] = 0.0;
      };
      return *this;
   }

   /*
   * Calculate component of vector v parallel to vector p.
   *
   * Upon return, the invoking vector is equal to the vector projection of
   * vector v along a direction parallel to vector p.
   */
   inline
   Vector& Vector::parallel(const Vector& v, const Vector& p)
   {
      double pp, fac;
      pp = p.square();
      if (pp > 1.0E-8) {
         fac = v.dot(p)/pp;
      } else {
         fac = 0.0;
      }
      elem_[0] = p.elem_[0]*fac;
      elem_[1] = p.elem_[1]*fac;
      elem_[2] = p.elem_[2]*fac;
      return *this;
   }

   /*
   * Calculate component of vector v transverse to vector p.
   *
   * Upon return, the invoking vector is equal to the vector
   * projection of vector v perpendicular to vector p.
   */
   inline
   Vector& Vector::transverse(const Vector& v, const Vector& p)
   {
      double pp, fac;
      pp = p.square();
      if (pp > 1.0E-8) {
         fac = v.dot(p)/pp;
      } else {
         fac = 0.0;
      }
      elem_[0] = v.elem_[0] - p.elem_[0]*fac;
      elem_[1] = v.elem_[1] - p.elem_[1]*fac;
      elem_[2] = v.elem_[2] - p.elem_[2]*fac;
      return *this;
   }

   /*
   * Compute the index corresponding to minimum element in a vector.
   */
   inline
   int minId(const Vector& v)
   {
      int id = 0;
      for (int i = 1; i < DimensionSq; i++) {
         if (v[i] < v[id])
            id = i;
      }
      return id;
   }

   /*
   * Compute the index corresponding to maximum element in a vector.
   */
   inline
   int maxId(const Vector& v)
   {
      int id = 0;
      for (int i = 1; i < DimensionSq; i++) {
         if (v[i] > v[id])
            id = i;
      }
      return id;
   }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   inline void Vector::serialize(Archive& ar, const unsigned int version)
   {
      ar & elem_[0];
      ar & elem_[1];
      ar & elem_[2];
   }

}
#endif
