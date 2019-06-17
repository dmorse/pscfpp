#ifndef PSCF_INT_VEC_H
#define PSCF_INT_VEC_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Vec.h"
#include <iostream>
#include <util/global.h>

namespace Pscf
{

   /**
   * An IntVec<D, T> is a D-component vector of elements of integer type T.
   *
   * Default of type T is T = int.
   *
   * \ingroup Pscf_Math_Module
   */
   template <int D, typename T = int>
   class IntVec : public Vec<D, T>
   {

   public:

      /// \name Constructors
      //@{

      /**
      * Default constructor
      */
      IntVec<D, T>()
        : Vec<D, T>()
      {}

      /**
      * Copy constructor
      *
      * \param v IntVec<D, T> to be copied
      */
      IntVec<D, T>(const IntVec<D, T>& v)
       : Vec<D, T>(v)
      {}

      /**
      * Construct from C array.
      *
      * \param v C array to be copied
      */
      IntVec<D, T>(T const * v)
       : Vec<D, T>(v)
      {}

      /**
      * Constructor, initialize all elements to a scalar value.
      *
      * \param s scalar initial value for all elements.
      */
      explicit IntVec<D, T>(T s)
       : Vec<D, T>(s)
      {}

      /// Width of field per Cartesian coordinate in stream IO
      static const int Width = 25;

   private:

      //template 
      // friend std::istream& operator >> (std::istream& in, IntVec<D, T> &vector);

      //template 
      //friend std::ostream& operator << (std::ostream& in, IntVec<D, T> &vector);

   };

   // Friend functions and operators
   
   /**
   * istream extractor for a IntVec<D, T>.
   *
   * Input elements of a vector from stream, without line breaks.
   *
   * \param in  input stream
   * \param vector  IntVec<D, T> to be read from stream
   * \return modified input stream
   */
   template <int D, typename T>
   std::istream& operator >> (std::istream& in, IntVec<D, T> &vector)
   {
      for (int i = 0; i < D; ++i) {
         in >> vector[i];
      }
      return in;
   }

   /**
   * ostream inserter for a IntVec<D, T>.
   *
   * Output a IntVec<D, T> to an ostream, without line breaks.
   *
   * Output elements of a vector to stream, without line breaks.
   * \param  out  output stream
   * \param  vector  IntVec<D, T> to be written to stream
   * \return modified output stream
   */
   template <int D, typename T>
   std::ostream& operator << (std::ostream& out, const IntVec<D, T> &vector) 
   {
      for (int i = 0; i < D; ++i) {
         out.width(IntVec<D, T>::Width);
         out << vector[i];
      }
      return out;
   }

   /**
   * Equality of two IntVec<D> objects.
   *
   * \return true if v1 == v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator == (const IntVec<D, T>& v1, const IntVec<D, T>& v2) 
   {
      for (int i = 0; i < D; ++i) {
         if (v1[i] != v2[i]) {
            return false;
         }
      }
      return true;
   }
   
   /**
   * Equality of an IntVec<D> and a Vec<D, T>
   *
   * \return true if v1 == v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator == (const IntVec<D, T>& v1, const Vec<D, T>& v2) 
   {
      for (int i = 0; i < D; ++i) {
         if (v1[i] != v2[i]) {
            return false;
         }
      }
      return true;
   }
   
   /**
   * Equality of an Vec<D, T> and an IntVec<D, T>
   *
   * \return true if v1 == v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator == (const Vec<D, T>& v1, const IntVec<D, T>& v2) 
   {  return (v2 == v1); }

   /**
   * Inequality of two IntVec<D, T> objects.
   *
   * \return true if v1 != v2, false if v1 == v2.
   */
   template <int D, typename T>
   inline
   bool operator != (const IntVec<D, T>& v1, const IntVec<D, T>& v2) 
   { return !(v1 == v2); }
   
   /**
   * Inequality of an IntVec<D> and a Vec<D, T>
   *
   * \return true if v1 == v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator != (const IntVec<D, T>& v1, const Vec<D, T>& v2) 
   {  return !(v1 == v2); }

   /**
   * Equality of an Vec<D, T> and an IntVec<D, T>
   *
   * \return true if v1 == v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator != (const Vec<D, T>& v1, const IntVec<D, T>& v2) 
   {  return !(v2 == v1); }

   /**
   * Less than comparison for two IntVec<D, T>s.
   *
   * \return true if v1 < v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator < (const IntVec<D, T>& v1, const IntVec<D, T>& v2) 
   {
      for (int i = 0; i < D; ++i) {
         if (v1[i] > v2[i]) {
           return false;
         } else 
         if (v1[i] < v2[i]) {
            return true;
         } 
      }
      return false;
   }
   
   /**
   * Greater than comparison for two IntVec<D, T>s.
   *
   * \return true if v1 > v2, false otherwise.
   */
   template <int D, typename T>
   inline 
   bool operator > (const IntVec<D, T>& v1, const IntVec<D, T>& v2) 
   {
      for (int i = 0; i < D; ++i) {
         if (v1[i] < v2[i]) {
           return false;
         } else 
         if (v1[i] > v2[i]) {
            return true;
         } 
      }
      return false;
   }
   
}
#endif
