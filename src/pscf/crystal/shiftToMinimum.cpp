/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "shiftToMinimum.h"

namespace Pscf
{

   using namespace Util;

   template <> 
   IntVec<1> shiftToMinimum(IntVec<1>& v, const IntVec<1> d, const UnitCell<1> cell)
   {
      IntVec<1> u;
      if (v[0] > d[0]/2) {
         u[0] = v[0] - d[0];
      } else {
         u[0] = v[0];
      }
      return u;
   }

   template <>
   IntVec<2> shiftToMinimum(IntVec<2>& v, const IntVec<2> d, const UnitCell<2> cell)
   {
      IntVec<2> u;
      for( int i = 0; i < 2; i++)
      {
         if (v[i] > d[i]/2) {
            u[i] = v[i] - d[i];
         } else {
            u[i] = v[i];
         }
      }
      return u;
   }

   template <>
   IntVec<3> shiftToMinimum(IntVec<3>& v, const IntVec<3> d, const UnitCell<3> cell)
   {
      IntVec<3> u;
      for( int i = 0; i < 3; i++)
      {
         if (v[i] > d[i]/2) {
            u[i] = v[i] - d[i];
         } else {
            u[i] = v[i];
         }
      }
      return u;
   }

}
