#ifndef UTIL_INT_VECTOR_H
#define UTIL_INT_VECTOR_H

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

namespace Util
{

   /**
   * An IntVector is an integer Cartesian vector.
   *
   * The Cartesian elements of a IntVector can be accessed using array notation:
   * The elements of a three dimensional IntVector v are v[0], v[1], and v[2].
   * The subscript operator [] returns elements as references, which can be
   * used on either the left or right side of an assignment operator.
   *
   * The arithmetic assignment operators +=, -=, *=, and /= are overloaded.
   * The operators += and -= represent increment or decrement by a vector,
   * while *= and /= represent multiplication or division by an integer.
   *
   * All other unary and binary mathematical operations are implemented as 
   * methods. Operations that yield a scalar result, such as a dot product, 
   * are implemented as methods that return the resulting value. Operations 
   * that yield a IntVector, such as vector addition, are implemented by methods 
   * that assign the result to the invoking vector, and return a reference 
   * to the invoking vector. For example,
   * \code
   *    IntVector a, b, c;
   *    int s;
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
   *    s = dot(a, b) 
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
   * This syntax for IntVector valued operations avoids dynamic allocation of
   * temporary IntVector objects, by requiring that the invoking function 
   * provide an object to hold the result.
   *
   * For efficiency, all methods in this class are inlined.
   *
   * \ingroup Space_Module
   */
   class IntVector
   {

   public:

      /// \name Constructors
      //@{ 
      
      /**
      * Default constructor
      */
      IntVector()
      {}

      /**
      * Copy constructor
      */
      IntVector(const IntVector& v)
      {
         elem_[0] = v.elem_[0];
         elem_[1] = v.elem_[1];
         elem_[2] = v.elem_[2];
      }

      /**
      * Constructor, initialize all elements to the same scalar.
      *
      * \param scalar initial value for all elements.
      */
      explicit IntVector(int scalar)
      {
         elem_[0] = scalar;
         elem_[1] = scalar;
         elem_[2] = scalar;
      }

      /**
      * Construct IntVector from C int[3] array.
      *
      * \param v array of 3 coordinates
      */
      explicit IntVector(const int* v)
      {
         elem_[0] = v[0];
         elem_[1] = v[1];
         elem_[2] = v[2];
      }

      /**
      * Construct IntVector from its coordinates.
      *
      * \param x x-axis coordinate
      * \param y y-axis coordinate
      * \param z z-axis coordinate
      */
      IntVector(int x, int y, int z = 0)
      {
         elem_[0] = x;
         elem_[1] = y;
         elem_[2] = z;
      }

      //@}
      
      /**
      * Set all elements of a 3D vector to zero.
      */
      IntVector& zero()
      {
         elem_[0] = 0;
         elem_[1] = 0;
         elem_[2] = 0;
         return *this;
      }

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

      /// \name Assignment
      //@{
      
      /**
      * Copy assignment.
      *
      * \param v IntVector to assign.
      */
      IntVector& operator=(const IntVector& v)
      {
         elem_[0] = v.elem_[0];
         elem_[1] = v.elem_[1];
         elem_[2] = v.elem_[2];
         return *this;
      }

      /**
      * Assignment from C int[] array.
      *
      * \param v array of coordinates
      */
      IntVector& operator=(const int* v)
      {
         elem_[0] = v[0];
         elem_[1] = v[1];
         elem_[2] = v[2];
         return *this;
      }

      //@}
      /// \name Arithmetic Assignment
      //@{
      
      /**
      * Add vector dv to this vector.
      *
      * Upon return, *this = this + dv.
      *
      * \param  dv   vector increment (input)
      */
      void operator+=(const IntVector& dv)
      {
         elem_[0] += dv.elem_[0];
         elem_[1] += dv.elem_[1];
         elem_[2] += dv.elem_[2];
      }

      /**
      * Subtract vector dv from this vector.
      *
      * Upon return, *this = this + dv.
      *
      * \param  dv   vector increment (input)
      */
      void operator-=(const IntVector& dv)
      {
         elem_[0] -= dv.elem_[0];
         elem_[1] -= dv.elem_[1];
         elem_[2] -= dv.elem_[2];
      }

      /**
      * Multiply this vector by scalar s.
      *
      * Upon return, *this = this*s.
      *
      * \param  s  scalar multiplier
      */
      void operator*=(int s)
      {
         elem_[0] *= s;
         elem_[1] *= s;
         elem_[2] *= s;
      }

      //@}
      /// \name Array Subscript 
      //@{
      
      /**
      * Return one Cartesian element by value.
      *
      * \param   i element index
      * \return  element i of the vector
      */
      const int& operator[](int i) const
      {
         assert(i >= 0);      
         assert(i < Dimension);      
         return elem_[i]; 
      }

      /**
      * Return a reference to one element of the vector.
      *
      * \param   i element index
      * \return  element i of the vector
      */
      int& operator[](int i)
      { 
         assert(i >= 0);      
         assert(i < Dimension);      
         return elem_[i]; 
      }

      //@}
      /// \name Scalar valued functions
      //@{
      
      /**
      * Return square magnitude of this vector.
      *
      * \return square magnitude of this vector
      */
      int square() const
      {  return elem_[0]*elem_[0] + elem_[1]*elem_[1] + elem_[2]*elem_[2]; }

      /**
      * Return dot product of this vector and vector v.
      *
      * \param v  input vector
      * \return dot product of this vector and vector v
      */
      int dot(const IntVector& v) const
      {
         return elem_[0]*v.elem_[0] + elem_[1]*v.elem_[1] + elem_[2]*v.elem_[2];
      }

      //@}
      /// \name IntVector valued functions (result assigned to invoking object)
      //@{

      /**
      * Add vectors v1 and v2.
      *
      * Upon return, *this = v1 + v2.
      *
      * \param v1  vector (input)
      * \param v2  vector (input)
      */
      IntVector& add(const IntVector& v1, const IntVector& v2)
      {
         elem_[0] = v1.elem_[0] + v2.elem_[0];
         elem_[1] = v1.elem_[1] + v2.elem_[1];
         elem_[2] = v1.elem_[2] + v2.elem_[2];
         return *this;
      }

      /**
      * Subtract vector v2 from v1.
      *
      * Upon return, *this = v1 - v2.
      *
      * \param v1  vector (input)
      * \param v2  vector (input)
      * \return modified invoking vector
      */
      IntVector& subtract(const IntVector& v1, const IntVector& v2)
      {
         elem_[0] = v1.elem_[0] - v2.elem_[0];
         elem_[1] = v1.elem_[1] - v2.elem_[1];
         elem_[2] = v1.elem_[2] - v2.elem_[2];
         return *this;
      }

      /**
      * Multiply a vector v by a scalar s.
      *
      * Upon return, *this = v*s.
      *
      * \param v  vector input
      * \param s  scalar input
      * \return modified invoking vector
      */
      IntVector& multiply(const IntVector& v, int s)
      {
         elem_[0] = v.elem_[0]*s;
         elem_[1] = v.elem_[1]*s;
         elem_[2] = v.elem_[2]*s;
         return *this;
      }

      /**
      * Calculate cross product of vectors v1 and v2.
      *
      * Upon return, *this = v1 x v2.
      *
      * \param v1  input  vector
      * \param v2  input  vector
      * \return modified invoking vector
      */
      IntVector& cross(const IntVector& v1, const IntVector& v2)
      {
         elem_[0] = v1.elem_[1]*v2.elem_[2] - v1.elem_[2]*v2.elem_[1];
         elem_[1] = v1.elem_[2]*v2.elem_[0] - v1.elem_[0]*v2.elem_[2];
         elem_[2] = v1.elem_[0]*v2.elem_[1] - v1.elem_[1]*v2.elem_[0];
         return *this;
      }

      //@}

      /// \name Static Members
      //@{ 
      
      /// Zero IntVector.
      static const IntVector Zero;
     
      /**
      * Initialize static IntVector::Zero.
      */
      static void initStatic();

      #ifdef UTIL_MPI
      /**
      * Commit MPI datatype MpiTraits<IntVector>::type.
      */
      static void commitMpiType();
      #endif

      //@}

   private:

      /// Width of field in stream IO
      static const int Width  = 15;

      /// Elements of the vector.
      int elem_[Dimension];

   //friends:

      friend bool operator == (const IntVector& v1, const IntVector& v2);

      friend bool operator == (const IntVector& v1, const int* v2);

      friend std::istream& operator >> (std::istream& in, IntVector &vector);

      friend std::ostream& operator << (std::ostream& out, const IntVector &vector);

   };

   /// Equality for IntVectors.
   bool operator == (const IntVector& v1, const IntVector& v2);

   /// Equality of IntVector and C array.
   bool operator == (const IntVector& v1, const int* v2);

   /// Equality of C array and IntVector.
   bool operator == (const int* v1, const IntVector& v2);

   /// Inequality of two IntVectors.
   bool operator != (const IntVector& v1, const IntVector& v2);

   /// Inequality of IntVector and C array.
   bool operator != (const IntVector& v1, const int* v2);

   /// Inequality of C array and IntVector.
   bool operator != (const int* v1, const IntVector& v2);

   /**
   * istream extractor for a IntVector.
   *
   * Input elements of a vector from stream, without line breaks.
   *
   * \param in      input stream
   * \param vector  IntVector to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, IntVector &vector);

   /**
   * ostream inserter for a IntVector.
   *
   * Output elements of a vector to stream, without line breaks.
   * \param  out     output stream
   * \param  vector  IntVector to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const IntVector &vector);

   #ifdef UTIL_MPI
   /**
   * Explicit specialization MpiTraits<IntVector>.
   */
   template <> class MpiTraits<IntVector> 
   {  
   public:  
      static MPI::Datatype type;    ///< MPI Datatype
      static bool hasType;          ///< Is the MPI type initialized?
   };
   #endif

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   inline 
   void IntVector::serialize(Archive& ar, const unsigned int version)
   { 
      ar & elem_[0];
      ar & elem_[1];
      ar & elem_[2];
   }

} 
#endif
