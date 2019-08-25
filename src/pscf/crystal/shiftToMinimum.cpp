/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "shiftToMinimum.h"
#include <util/global.h>

namespace Pscf
{

   using namespace Util;

   template <> 
   IntVec<1> shiftToMinimum(IntVec<1>& v, const IntVec<1> d, const UnitCell<1> cell)
   {
      IntVec<1> u;
      if (v[0] > d[0]/2) {
         u[0] = v[0] - d[0];
      } else 
      if (v[0] < -(d[0]/2) ) {
         u[0] = v[0] + d[0];
      } else {
         u[0] = v[0];
      }
      UTIL_ASSERT(abs(u[0]) <= d[0]/2);
      return u;
   }

   template <>
   IntVec<2> shiftToMinimum(IntVec<2>& v, const IntVec<2> d, 
                            const UnitCell<2> cell)
   {
      double Gsq;
      double Gsq_min = cell.ksq(v);
      IntVec<2> r = v;               // Minimum image
      IntVec<2> u;                   // Vector used in loop
      int s0, s1;
      for (s0 = 1; s0 >= - 1; --s0) { 
         u[0] = v[0] + s0*d[0];
         for (s1 = 1; s1 >= - 1; --s1) {
            u[1] = v[1] + s1*d[1];
            Gsq = cell.ksq(u);
            if (Gsq < Gsq_min) {
               Gsq_min = Gsq;
               r = u;
            }
         }
      }
      return r;
   }

   template <>
   IntVec<3> shiftToMinimum(IntVec<3>& v, const IntVec<3> d, 
                            const UnitCell<3> cell)
   {
      double Gsq;
      double Gsq_min = cell.ksq(v);
      IntVec<3> r = v;               // Minimum image
      IntVec<3> u;                   // Vector used in loop
      int s0, s1, s2;
      for (s0 = 1; s0 >= -1; --s0) { 
         u[0] = v[0] + s0*d[0];
         for (s1 = 1; s1 >= -1; --s1) {
            u[1] = v[1] + s1*d[1];
            for (s2 = 1; s2 >= -1; --s2) {
               u[2] = v[2] + s2*d[2];
               Gsq = cell.ksq(u);
               if (Gsq < Gsq_min) {
                  Gsq_min = Gsq;
                  r = u;
               }
            }
         }
      }
      return r;
   }

}
