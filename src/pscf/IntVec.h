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

   // Forward declarations

   template <int D, typename T> class IntVec;

   template <int D, typename T>
   bool operator == (const IntVec<D, T>& v1, const IntVec<D, T>& v2);

   template <int D, typename T>
   bool operator == (const IntVec<D, T>& v1, const T* v2);

   template <int D, typename T>
   std::istream& operator >> (std::istream& in, IntVec<D, T> &vector);

   template <int D, typename T>
   std::ostream& operator << (std::ostream& out, const IntVec<D, T> &vector);

   /**
   * A IntVec<D, T> is D-component vector with elements of integer type.
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
      IntVec<D, T>();

      /**
      * Copy constructor
      *
      * \param v IntVec<D, T> to be copied
      */
      IntVec<D, T>(const IntVec<D, T>& v);

      /**
      * Constructor, initialize all elements to a scalar value.
      *
      * \param scalar initial value for all elements.
      */
      explicit IntVec<D, T>(T scalar);

      //@}
      #if 0
      /// \name Static Members
      //@{

      /**
      * Zero IntVec<D, T> = {0.0, 0.0, 0.0}
      */
      static const IntVec<D, T> Zero;

      //@}
      #endif

   private:

      /// Width of field per Cartesian coordinate in stream IO
      static const int Width = 25;

      /// Precision in stream IO of IntVec<D, T> coordinates
      static const int Precision = 17;

   //friends:

      friend 
      bool operator == <>(const IntVec<D, T>& v1, const IntVec<D, T>& v2);

      friend 
      bool operator == <>(const IntVec<D, T>& v1, const T* v2);

      friend std::istream& 
      operator >> <>(std::istream& in, IntVec<D, T> &vector);

      friend std::ostream& 
      operator << <>(std::ostream& out, const IntVec<D, T> &vector);

   };

   /// Equality for IntVec<D, T>s.
   template <int D, typename T> 
   bool operator == (const IntVec<D, T>& v1, const IntVec<D, T>& v2);

   /// Inequality of two IntVec<D, T>s.
   template <int D, typename T> 
   bool operator != (const IntVec<D, T>& v1, const IntVec<D, T>& v2);

   /**
   * istream extractor for a IntVec<D, T>.
   *
   * Input elements of a vector from stream, without line breaks.
   *
   * \param in      input stream
   * \param vector  IntVec<D, T> to be read from stream
   * \return modified input stream
   */
   template <int D, typename T> 
   std::istream& operator >> (std::istream& in, IntVec<D, T> &vector);

   /**
   * ostream inserter for a IntVec<D, T>.
   *
   * Output elements of a vector to stream, without line breaks.
   * \param  out     output stream
   * \param  vector  IntVec<D, T> to be written to stream
   * \return modified output stream
   */
   template <int D, typename T> 
   std::ostream& operator << (std::ostream& out, const IntVec<D, T> &vector);

   // Inline methods

   /*
   * Default constructor
   */
   template <int D, typename T> 
   inline IntVec<D, T>::IntVec()
    : Vec<D,T>()
   {}

   /*
   * Copy constructor
   */
   template <int D, typename T> 
   inline IntVec<D, T>::IntVec(const IntVec<D, T>& v)
    : Vec<D,T>(v)
   {}

   /*
   * Constructor, initialize all elements to a scalar value.
   */
   template <int D, typename T> 
   inline IntVec<D, T>::IntVec(T s)
    : Vec<D,T>(s)
   {}

}
#endif
