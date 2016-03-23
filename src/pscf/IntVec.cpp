/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntVec.h"
#include <util/global.h>

namespace Pscf
{

   #if 0
   // Define static Zero vector
   const IntVec<D, T> IntVec<D, T>::Zero = IntVec<D, T>(0);
   #endif

   // Equality and inequality operators
   
   template <int D, typename T>
   bool operator==(const IntVec<D, T>& v1, const IntVec<D, T>& v2) 
   {
      for (int i=0; i < D; ++i) {
         if (v1.elem_[i] != v2.elem_[i]) {
            return false;
         }
      }
      return true;
   }
   
   template <int D, typename T>
   bool operator!=(const IntVec<D, T>& v1, const IntVec<D, T>& v2) 
   { return !(v1 == v2); }
   
   /* 
   * Input a IntVec<D, T> from an istream, without line breaks.
   */
   template <int D, typename T>
   std::istream& operator>>(std::istream& in, IntVec<D, T> &vector)
   {
      for (int i=0; i < D; ++i) {
         in >> vector.elem_[i];
      }
      return in;
   }
   
   /* 
   * Output a IntVec<D, T> to an ostream, without line breaks.
   */
   template <int D, typename T>
   std::ostream& operator<<(std::ostream& out, const IntVec<D, T> &vector) 
   {
      for (int i=0; i < D; ++i) {
         out.width(IntVec<D, T>::Width);
         out << vector.elem_[i];
      }
      return out;
   }

} 
