#ifndef PSCF_VEC_H
#define PSCF_VEC_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/accumulators/setToZero.h>
#include <iostream>

using namespace Util;

namespace Pscf
{

   /**
   * A Vec<D, T><D,T> is a D-component vector with elements of type T.
   *
   * The elements of a Vec<D, T> can be accessed using subscript operator, 
   * as for a built in array.
   *
   * The arithmetic assignment operators +=, -=, and *= are overloaded 
   * to allow vector-vector addition and subtraction and vector-scalar
   * multiplication.
   *
   * All other unary and binary mathematical operations are implemented 
   * as methods or free functions. Operations that yield a Vec<D, T>, such 
   * as addition, are implemented by methods that assign the result to the 
   * invoking Vec object, and return this object by reference. For example,
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
   * This syntax for functions that yield a vector makes the allocation of 
   * a temporary Vec<D, T> object explicit, by requiring that the invoking 
   * function be a member of an object that will hold the result.
   *
   * For efficiency, all member functions are declared inline.
   *
   * \ingroup Pscf_Math_Module
   */
   template <int D, typename T>
   class Vec
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
      * Constructor from a C-array.
      *
      * \param v array to be copied
      */
      explicit Vec<D, T>(T const *  v);

      /**
      * Constructor, initialize all elements to a common scalar value.
      *
      * \param s initial value for all elements.
      */
      explicit Vec<D, T>(T s);

      //@}

      /// \name Assignment and Initialization
      //@{

      /**
      * Copy assignment.
      *
      * \param v Vec<D, T> to assign.
      * \return this object, after modification
      */
      Vec<D, T>& operator = (const Vec<D, T>& v);

      /**
      * Assignment all elements to the same scalar T value.
      *
      * \param s scalar value
      * \return this object, after modification
      */
      Vec<D, T>& operator = (T s);

      /**
      * Set all elements to zero.
      *
      * \return this object, after modification
      */
      Vec<D, T>& setToZero();

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
      * Add a common scalar to all components.
      *
      * \param s  scalar additive constant (input)
      */
      void operator += (T s);

      /**
      * Subtract a common scalar from all components.
      *
      * \param s  scalar subtractive constant (input)
      */
      void operator -= (T s);

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
      /// \name Vec<D, T> valued functions (assigned to invoking object)
      //@{

      /**
      * Add vectors v1 and v2.
      *
      * Upon return, *this = v1 + v2.
      *
      * \param v1 vector (input)
      * \param v2 vector (input)
      * \return modified invoking vector
      */
      Vec<D, T>& add(const Vec<D, T>& v1, const Vec<D, T>& v2);

      /**
      * Subtract vector v2 from v1.
      *
      * Upon return, *this = v1 - v2.
      *
      * \param v1 vector (input)
      * \param v2 vector (input)
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

      /**
      * Return negative of vector v.
      *
      * Upon return, *this = -v;
      *
      * \param v  vector input
      * \return modified invoking vector
      */
      Vec<D, T>& negate(const Vec<D, T>& v);

      /**
      * Negate all elements of this vector.
      *
      * Upon return, all elements of this have been negated (reversed)
      *
      * \return this object, after modification
      */
      Vec<D, T>& negate();

      //@}

      /**
      * Serialize to/from an archive.
      *
      * Implementation uses syntax of Boost::serialize.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Width of field per Cartesian coordinate in stream IO
      static const int Width = 25;

      /// Precision in stream IO of Vec<D, T> coordinates
      static const int Precision = 17;

      /// Elements of the vector.
      T elem_[D];

   };

   // Associated functions

   /**
   * Return dot product of two vectors.
   *
   * \param v1 first input vector
   * \param v2 second input vector
   * \return dot product v1.v2
   */
   template <int D, typename T>
   inline 
   T dot(Vec<D, T> const & v1, Vec<D, T> const & v2)
   {
      T value;
      setToZero(value);
      for (int i = 0; i < D; ++i) {
         value += v1[i]*v2[i];
      }
      return value;
   }

   /**
   * Return the sum of two vectors.
   *
   * \param v1 first input vector
   * \param v2 second input vector
   * \return sum v1 + v2
   */
   template <int D, typename T>
   inline 
   Vec<D, T> operator + (Vec<D, T> const & v1, Vec<D, T> const & v2)
   {
      Vec<D, T> value;
      value.add(v1, v2);
      return value;
   }

   // Inline method definitions

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
   * Constructor, from C-Array.
   */
   template <int D, typename T> 
   inline Vec<D, T>::Vec(T const * v)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = v[i];
      }
   }

   /*
   * Constructor, initialize all elements to a scalar value s.
   */
   template <int D, typename T> 
   inline Vec<D, T>::Vec(T s)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = s;
      }
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
   * Assign all elements to a common scalar value.
   */
   template <int D, typename T> 
   inline Vec<D, T>& Vec<D, T>::operator = (T s)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = s;
      }
      return *this;
   }

   /*
   * Set all elements of a 3D vector to zero.
   */
   template <int D, typename T> 
   inline Vec<D, T>& Vec<D, T>::setToZero()
   {
      for (int i = 0; i < D; ++i) {
         setToZero(elem_[i]); 
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
   * Add a common scalar to all components.
   */
   template <int D, typename T> 
   inline void Vec<D, T>::operator += (T s)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] += s;
      }
   }

   /*
   * Subtract a common scalar from all components.
   */
   template <int D, typename T> 
   inline void Vec<D, T>::operator -= (T s)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] -= s;
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
      assert(i < D); 
      return elem_[i]; 
   }

   /*
   * Return a reference to one element of the vector.
   */
   template <int D, typename T>
   inline T& Vec<D, T>::operator [] (int i)
   {
      assert(i >=0); 
      assert(i < D); 
      return elem_[i]; 
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
   * Compute and return negation of a vector.
   *
   * Upon return, *this = -v.
   */
   template <int D, typename T>
   inline 
   Vec<D, T>& Vec<D, T>::negate(Vec<D, T> const & v)
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = -v.elem_[i];
      }
      return *this;
   }

   /*
   * Negate (reverse sign) of this vector.
   *
   * Upon return, *this = -v.
   */
   template <int D, typename T>
   inline 
   Vec<D, T>& Vec<D, T>::negate()
   {
      for (int i = 0; i < D; ++i) {
         elem_[i] = -elem_[i];
      }
      return *this;
   }

   /*
   * Serialize to/from an archive.
   */
   template <int D, typename T>
   template <class Archive>
   inline 
   void Vec<D, T>::serialize(Archive& ar, const unsigned int version)
   { 
      for (int i = 0; i < D; ++i) {
         ar & elem_[i];
      }
   }

}
#endif
