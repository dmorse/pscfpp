/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntVector.h"
#include "Dimension.h"
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

namespace Util
{

   // Define static Zero vector
   const IntVector IntVector::Zero = IntVector(0);

   // Equality operators
   
   bool operator==(const IntVector& v1, const IntVector& v2) 
   {
      for (int i=0; i < Dimension; ++i) {
         if (v1.elem_[i] != v2.elem_[i]) {
            return false;
         }
      }
      return true;
   }
   
   
   bool operator==(const IntVector& v1, const int* v2) 
   {
      for (int i=0; i < Dimension; ++i) {
         if (v1.elem_[i] != v2[i]) {
            return false;
         }
      }
      return true;
   }
   
   bool operator==(const int* v1, const IntVector& v2) 
   { return (v2 == v1); }
   
   // Inequality operators
   
   bool operator!=(const IntVector& v1, const IntVector& v2) 
   { return !(v1 == v2); }
   
   
   bool operator!=(const IntVector& v1, const int* v2) 
   { return !(v1 == v2); }
   
   
   bool operator!=(const int* v1, const IntVector& v2) 
   { return !(v2 == v1); }
   
   /* 
   * Input a IntVector from an istream, without line breaks.
   */
   std::istream& operator>>(std::istream& in, IntVector &vector)
   {
      for (int i=0; i < Dimension; ++i) {
         in >> vector.elem_[i];
      }
      return in;
   }
   
   /* 
   * Output a IntVector to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const IntVector &vector) 
   {
      for (int i=0; i < Dimension; ++i) {
         out.width(IntVector::Width);
         out << vector.elem_[i];
      }
      return out;
   }

   #ifdef UTIL_MPI
   // Initialize MpiTraits<IntVector>
   MPI::Datatype MpiTraits<IntVector>::type = MPI::BYTE;
   bool MpiTraits<IntVector>::hasType = false;

   /*
   * Commit MPI Datatype.
   */
   void IntVector::commitMpiType() 
   {
      if (!MpiTraits<IntVector>::hasType) {
         MpiStructBuilder builder;
         IntVector vector;
         builder.setBase(&vector);
         builder.addMember(&vector[0], MPI::INT);
         builder.addMember(&vector[1], MPI::INT);
         builder.addMember(&vector[2], MPI::INT);
         builder.commit(MpiTraits<IntVector>::type);
         MpiTraits<IntVector>::hasType = true;
      }
   }
   #endif

   /*
   * This static method exists to guarantee initialization of static 
   * constants and variables defined in this file.  Call it somewhere 
   * in the program to guarantee that the contents of this file will 
   * be linked, rather than optimized away. It may only be called once.
   */
   void IntVector::initStatic() 
   {
      static int nCall = 0;
      ++nCall;
   }

} 
