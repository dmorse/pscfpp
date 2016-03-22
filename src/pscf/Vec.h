#ifndef PSCF_VEC_H
#define PSCF_VEC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <iostream>
#include <cmath>

namespace Pscf
{

   /**
   * A Vec<D, T><D,T> is a simple d-component vector with elements of type T.
   *
   * The elements of a Vec<D, T> can be accessed using subscript operator,
   * as for a built in array.
   *
   * The arithmetic assignment operators +=, -=, and *= are overloaded.
   *
   * All other unary and binary mathematical operations are implemented as
   * methods or friend functions. Operations that yield a Vec<D, T>, such 
   * as addition, are implemented by methods that assign the result to the 
   * invoking vector, and return a reference to the invoking vector. For 
   * example,
   * \code
   *
   *    Vec<3, double> a, b, c;
   *
   *    a[0] = 0.0
   *    a[1] = 1.0
   *    a[2] = 2.0
   *
   *    b[0] =  0.5
   *    b[1] = -0.5
   *    b[2] = -1.5
   *
   *    // Set a = a + b
   *    a += b
   *
   *    // Set b = b*2
   *    b *= 2.0;
   *
   *    // Set c = a + b
   *    c.add(a, b);
   *
   * \endcode
   * This syntax for Vec<D, T> valued operations avoids dynamic allocation 
   * of temporary Vec<D, T> objects, by requiring that the invoking 
   * function provide an object to hold the result.
   *
   * For efficiency, all methods in this class are declared inline.
   */
   template <int D, typename T>
   class Vec<D, T>
   {

   public:

      /// \name Constructors
      //@{

      /**
      * Default constructor
      */
      Vec<D, T>();

      /**
      * Copy constructor
      *
      * \param v Vec<D, T> to be copied
      */
      Vec<D, T>(const Vec<D, T>& v);

      /**
      * Constructor, initialize all elements to a scalar value.
      *
      * \param scalar initial value for all elements.
      */
      explicit Vec<D, T>(T scalar);

      //@}

      /**
      * Set all elements to zero.
      */
      Vec<D, T>& setToZero();

      /// \name Assignment
      //@{

      /**
      * Copy assignment.
      *
      * \param v Vec<D, T> to assign.
      */
      Vec<D, T>& operator = (const Vec<D, T>& v);

      /**
      * Assignment from C T[3] array.
      *
      * \param v  C-array of components
      */
      Vec<D, T>& operator = (const T* v);

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
      void operator += (const Vec<D, T>& dv);

      /**
      * Subtract vector dv from this vector.
      *
      * Upon return, *this = this + dv.
      *
      * \param dv  vector increment (input)
      */
      void operator -= (const Vec<D, T>& dv);

      /**
      * Multiply this vector by scalar s.
      *
      * Upon return, *this = (*this)*s.
      *
      * \param s  scalar multiplier
      */
      void operator *= (T s);

      //@}
      /// \name Array Subscript
      //@{

      /**
      * Return one Cartesian element by value.
      *
      * \param i  element index
      * \return element i of the vector
      */
      const T& operator [] (int i) const;

      /**
      * Return one element of the vector by references.
      *
      * \param i  element index
      * \return element i of this vector
      */
      T& operator [] (int i);

      //@}
      /// \name Scalar-valued functions
      //@{

      /**
      * Return dot product of this vector and vector v.
      *
      * \param  v input vector
      * \return dot product of this vector and vector v
      */
      T dot(const Vec<D, T>& v) const;

      //@}
      /// \name Vec<D, T> valued functions (assigned to invoking object)
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
      Vec<D, T>& add(const Vec<D, T>& v1, const Vec<D, T>& v2);

      /**
      * Subtract vector v2 from v1.
      *
      * Upon return, *this = v1 - v2.
      *
      * \param v1  vector (input)
      * \param v2  vector (input)
      * \return modified invoking vector
      */
      Vec<D, T>& subtract(const Vec<D, T>& v1, const Vec<D, T>& v2);

      /**
      * Multiply a vector v by a scalar s.
      *
      * Upon return, *this = v*s.
      *
      * \param v  vector input
      * \param s  scalar input
      * \return modified invoking vector
      */
      Vec<D, T>& multiply(const Vec<D, T>& v, T s);

      //@}
      /// \name Static Members
      //@{

      /**
      * Zero Vec<D, T> = {0.0, 0.0, 0.0}
      */
      static const Vec<D, T> Zero;

      //@}

   private:

      /// Width of field per Cartesian coordinate in stream IO
      static const int Width = 25;

      /// Precision in stream IO of Vec<D, T> coordinates
      static const int Precision = 17;

      /// Elements of the vector.
      T elem_[Dimension];

   //friends:

      friend bool operator == (const Vec<D, T>& v1, const Vec<D, T>& v2);

      friend bool operator == (const Vec<D, T>& v1, const T* v2);

      friend 
      std::istream& operator >> (std::istream& in, Vec<D, T> &vector);

      friend 
      std::ostream& operator << (std::ostream& out, const Vec<D, T> &vector);

   };

   /// Equality for Vec<D, T>s.
   bool operator == (const Vec<D, T>& v1, const Vec<D, T>& v2);

   /// Equality of Vec<D, T> and C array.
   bool operator == (const Vec<D, T>& v1, const T* v2);

   /// Inequality of two Vec<D, T>s.
   bool operator != (const Vec<D, T>& v1, const Vec<D, T>& v2);

   /// Inequality of Vec<D, T> and C array.
   bool operator != (const Vec<D, T>& v1, const T* v2);

   /// Inequality of C array and Vec<D, T>.
   bool operator != (const T* v1, const Vec<D, T>& v2);

   /**
   * istream extractor for a Vec<D, T>.
   *
   * Input elements of a vector from stream, without line breaks.
   *
   * \param in      input stream
   * \param vector  Vec<D, T> to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Vec<D, T> &vector);

   /**
   * ostream inserter for a Vec<D, T>.
   *
   * Output elements of a vector to stream, without line breaks.
   * \param  out     output stream
   * \param  vector  Vec<D, T> to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Vec<D, T> &vector);

   // Inline methods

   /*
   * Default constructor
   */
   template <int D, typename T> 
   inline Vec<D, T>::Vec()
   {}

   /*
   * Copy constructor
   */
   template <int D, typename T> 
   inline Vec<D, T>::Vec(const Vec<D, T>& v)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v.elem_[i];
      }
   }

   /*
   * Constructor, initialize all elements to a scalar value.
   */
   template <int D, typename T> 
   inline Vec<D, T>::Vec(T scalar)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = scalar;
      }
   }

   /*
   * Set all elements of a 3D vector to zero.
   */
   template <int D, typename T> 
   inline Vec<D, T>& Vec<D, T>::zero()
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = 0;
      }
      return *this;
   }

   /*
   * Assignment.
   */
   template <int D, typename T> 
   inline Vec<D, T>& Vec<D, T>::operator = (const Vec<D, T>& v)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v.elem_[i];
      }
      return *this;
   }

   /*
   * Assignment from C T[] array.
   */
   template <int D, typename T> 
   inline Vec<D, T>& Vec<D, T>::operator = (const T* v)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v[i];
      }
      return *this;
   }

   /*
   * Add vector dv to this vector.
   */
   template <int D, typename T> 
   inline void Vec<D, T>::operator += (const Vec<D, T>& dv)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] += dv.elem_[i];
      }
   }

   /*
   * Subtract vector dv from this vector.
   */
   template <int D, typename T> 
   inline void Vec<D, T>::operator -= (const Vec<D, T>& dv)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] -= dv.elem_[i];
      }
   }

   /*
   * Multiply this vector by scalar s.
   */
   template <int D, typename T> 
   inline void Vec<D, T>::operator *= (T s)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] *= s;
      }
   }

   /*
   * Return one Cartesian element by value.
   */
   template <int D, typename T>
   inline const T& Vec<D, T>::operator [] (int i) const
   { 
      assert(i >=0); 
      assert(i < Dimension); 
      return elem_[i]; 
   }

   /*
   * Return a reference to one element of the vector.
   */
   */
   template <int D, typename T>
   inline T& Vec<D, T>::operator [] (int i)
   {
      assert(i >=0); 
      assert(i < Dimension); 
      return elem_[i]; 
   }

   /*
   * Return square magnitude of this vector.
   */
   inline
   T Vec<D, T>::square() const
   {
      return (elem_[0]*elem_[0] + elem_[1]*elem_[1] + elem_[2]*elem_[2]);
   }

   /*
   * Return absolute magnitude of this vector.
   */
   template <int D, typename T>
   inline T Vec<D, T>::abs() const
   {  return sqrt(square()); }

   /*
   * Return dot product of this vector and vector v.
   */
   template <int D, typename T>
   inline T Vec<D, T>::dot(Vec<D, T> const & v) const
   {
      T value = 0;
      for (i = 0; i < D; ++i) {
         value += elem_[i]*v.elem_[i];
      }
      return value;
   }

   /*
   * Add vectors v1 and v2.
   *
   * Upon return, *this = v1 + v2.
   */
   template <int D, typename T>
   inline
   Vec<D, T>& Vec<D, T>::add(Vec<D, T> const & v1, 
                             Vec<D, T> const & v2)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v1.elem_[i] + v2.elem_[i];
      }
      return *this;
   }

   /*
   * Subtract vector v2 from v1.
   *
   * Upon return, *this = v1 - v2.
   * \return modified invoking vector
   */
   template <int D, typename T>
   inline
   Vec<D, T>& Vec<D, T>::subtract(Vec<D, T> const & v1, Vec<D, T> const & v2)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v1.elem_[i] - v2.elem_[i];
      }
      return *this;
   }

   /*
   * Multiply a vector v by a scalar s.
   *
   * Upon return, *this = v*s.
   */
   template <int D, typename T>
   inline 
   Vec<D, T>& Vec<D, T>::multiply(Vec<D, T> const & v, T s)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v.elem_[i]*s;
      }
      return *this;
   }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   inline void Vec<D, T>::serialize(Archive& ar, const unsigned int version)
   {
      for (int i=0; i < D; ++i) {
         ar & elem_[i];
      }
   }

}
#endif
