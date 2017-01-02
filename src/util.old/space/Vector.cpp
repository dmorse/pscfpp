/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Vector.h"
#include "Dimension.h"
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

namespace Util
{

   // Define the Zero Vector
   const Vector Vector::Zero = Vector(0.0);

   // Equality operators
   
   #define UTIL_VECTOR_EPSILON 1.0E-8
   
   bool operator==(const Vector& v1, const Vector& v2) 
   {
      for (int i=0; i < Dimension; ++i) {
         if ( fabs(v1.elem_[i] - v2.elem_[i]) > UTIL_VECTOR_EPSILON) {
            return false;
         }
      }
      return true;
   }
   
   bool operator==(const Vector& v1, const double* v2) 
   {
      for (int i=0; i < Dimension; ++i) {
         if ( fabs(v1.elem_[i] - v2[i]) > UTIL_VECTOR_EPSILON) {
            return false;
         }
      }
      return true;
   }
   
   #undef UTIL_VECTOR_EPSILON
   
   bool operator==(const double* v1, const Vector& v2) 
   { return (v2 == v1); }
   
   // Inequality operators
   
   bool operator!=(const Vector& v1, const Vector& v2) 
   { return !(v1 == v2); }
   
   bool operator!=(const Vector& v1, const double* v2) 
   { return !(v1 == v2); }
   
   bool operator!=(const double* v1, const Vector& v2) 
   { return !(v2 == v1); }
   
   /* 
   * Extract a Vector from an istream.
   */
   std::istream& operator>>(std::istream& in, Vector &vector)
   {
      for (int i=0; i < Dimension; ++i) {
         in >> vector.elem_[i];
      }
      return in;
   }
   
   /* 
   * Output a Vector to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const Vector &vector) 
   {
      for (int i=0; i < Dimension; ++i) {
         out.setf(std::ios::scientific);
         out.width(Vector::Width);
         out.precision(Vector::Precision);
         out << vector.elem_[i];
      }
      return out;
   }

   #ifdef UTIL_MPI
   /**
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<Vector>::type = MPI::BYTE;
   bool MpiTraits<Vector>::hasType = false;

   /**
   * Commit MPI Datatype.
   */
   void Vector::commitMpiType() 
   {
      if (!MpiTraits<Vector>::hasType) {
         MpiStructBuilder builder;
         Vector vector;
   
         builder.setBase(&vector);
         builder.addMember(&vector[0], MPI::DOUBLE);
         builder.addMember(&vector[1], MPI::DOUBLE);
         builder.addMember(&vector[2], MPI::DOUBLE);
         builder.commit(MpiTraits<Vector>::type);
         MpiTraits<Vector>::hasType = true;
      }
   }
   #endif

   /*
   * This static method exists to guarantee initialization of the static 
   * constant Vector::Zero that is defined in this file.  Call it once
   * in the program to guarantee that the contents of this file will be
   * linked, rather than optimized away. 
   */
   void Vector::initStatic()
   {
      static int nCall = 0;
      ++nCall;
   }

} 
