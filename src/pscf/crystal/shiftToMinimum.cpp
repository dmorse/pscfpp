/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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
      // Initialize minimum to input value
      const double epsilon = 1.0E-6;
      double Gsq = cell.ksq(v);
      double Gsq_min_lo = Gsq - epsilon;
      double Gsq_min_hi = Gsq + epsilon;

      // Loop over neighboring images 
      IntVec<2> r = v;                    // Minimum image of v
      IntVec<2> u;                        // Vector used in loop
      int s0, s1;
      const int q = 2;
      for (s0 = -q; s0 < q+1; ++s0) { 
         u[0] = v[0] + s0*d[0];
         for (s1 = -q; s1 < q+1; ++s1) {
            u[1] = v[1] + s1*d[1];
            Gsq = cell.ksq(u);
            if (Gsq < Gsq_min_hi) {
               if ((Gsq < Gsq_min_lo) || (r < u)) {
                  Gsq_min_lo = Gsq - epsilon;
                  Gsq_min_hi = Gsq + epsilon;
                  r = u;
               }
            }
         }
      }
      return r;
   }

   template <>
   IntVec<3> shiftToMinimum(IntVec<3>& v, const IntVec<3> d, 
                            const UnitCell<3> cell)
   {
      // Initialize minimum to input value
      const double epsilon = 1.0E-6;
      double Gsq = cell.ksq(v);
      double Gsq_min_lo = Gsq - epsilon;
      double Gsq_min_hi = Gsq + epsilon;

      // Loop over neighboring images
      IntVec<3> r = v;               // Minimum image
      IntVec<3> u;                   // Vector used in loop
      int s0, s1, s2;
      const int q = 2;
      for (s0 = -q; s0 < q + 1; ++s0) { 
         u[0] = v[0] + s0*d[0];
         for (s1 = -q; s1 < q + 1; ++s1) {
            u[1] = v[1] + s1*d[1];
            for (s2 = -q; s2 < q + 1; ++s2) {
               u[2] = v[2] + s2*d[2];
               Gsq = cell.ksq(u);
               if (Gsq < Gsq_min_hi) {
                  if ((Gsq < Gsq_min_lo) || (r < u)) {
                     Gsq_min_lo = Gsq - epsilon;
                     Gsq_min_hi = Gsq + epsilon;
                     r = u;
                  }
               }
            }
         }
      }
      return r;
   }

}
