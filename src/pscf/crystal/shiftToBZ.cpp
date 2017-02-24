/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "shiftToBZ.h"

namespace Pscf
{

   using namespace Util;

   template <> 
   IntVec<1> UnitCell(IntVec<1>& v, IntVec<1> d, UnitCell<1> cell)
   {
      IntVec<1> u;
      if (v[0] > d[0]/2) {
         u[0] = v[0] - d[0];
      } else {
         u[0] = v[0];
      }
      return u;
   }

}
