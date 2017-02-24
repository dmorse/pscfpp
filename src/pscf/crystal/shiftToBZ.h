#ifndef PSCF_SHIFT_BZ_H
#define PSCF_SHIFT_BZ_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <iostream>

namespace Pscf
{

   using namespace Util;

   template <int D> 
   IntVec<D> UnitCell(IntVec<D>& v, IntVec<D> d, UnitCell<D> cell);

   template <> 
   IntVec<1> UnitCell(IntVec<1>& v, IntVec<1> d, UnitCell<1> cell)
   #if 0
   {
      IntVec<1> u;
      if (v[0] > d[0]/2) {
         u[0] = v[0] - d[0];
      } else {
         u[0] = v[0];
      }
      return u;
   }
   #endif

}
#endif
