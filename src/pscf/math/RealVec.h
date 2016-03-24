#ifndef PSCF_REAL_VEC_H
#define PSCF_REAL_VEC_H

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
   * A RealVec<D, T> is D-component vector with elements of floating type T.
   */
   template <int D, typename T = double>
   class RealVec : public Vec<D, T>
   {

   public:

      /// \name Constructors
      //@{

      /**
      * Default constructor
      */
      RealVec<D, T>()
        : Vec<D, T>()
      {}

      /**
      * Copy constructor
      *
      * \param v RealVec<D, T> to be copied
      */
      RealVec<D, T>(const RealVec<D, T>& v)
       : Vec<D, T>(v)
      {}

      /**
      * Construct from C array.
      *
      * \param v C array to be copied
      */
      RealVec<D, T>(T const * v)
       : Vec<D, T>(v)
      {}

      /**
      * Constructor, initialize all elements to a scalar value.
      *
      * \param s scalar initial value for all elements.
      */
      explicit RealVec<D, T>(T s)
       : Vec<D, T>(s)
      {}

   private:

      /// Width of field per Cartesian coordinate in stream IO
      static const int Width = 25;

      /// Precision in stream IO of RealVec<D, T> coordinates
      static const int Precision = 17;

   };

   // Friend functions and operators

   /**
   * istream extractor for a RealVec<D, T>.
   *
   * Input elements of a vector from stream, without line breaks.
   *
   * \param in  input stream
   * \param vector  RealVec<D, T> to be read from stream
   * \return modified input stream
   */
   template <int D, typename T>
   std::istream& operator>>(std::istream& in, RealVec<D, T> &vector)
   {
      for (int i = 0; i < D; ++i) {
         in >> vector[i];
      }
      return in;
   }

   /**
   * ostream inserter for a RealVec<D, T>.
   *
   * Output a RealVec<D, T> to an ostream, without line breaks.
   *
   * Output elements of a vector to stream, without line breaks.
   * \param  out  output stream
   * \param  vector  RealVec<D, T> to be written to stream
   * \return modified output stream
   */
   template <int D, typename T>
   std::ostream& operator<<(std::ostream& out, const RealVec<D, T> &vector)
   {
      for (int i = 0; i < D; ++i) {
         out.width(RealVec<D, T>::Width);
         out << vector[i];
      }
      return out;
   }

}
#endif
