/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "shiftToMinimum.h"

namespace Pscf{
namespace Pssp_gpu{

   using namespace Util;

   template <> 
   IntVec<1> shiftToMinimum(IntVec<1>& v, const IntVec<1> d);

   template <> 
   IntVec<2> shiftToMinimum(IntVec<2>& v, const IntVec<2> d);

   template <> 
   IntVec<3> shiftToMinimum(IntVec<3>& v, const IntVec<3> d);

   template<>
   void shiftToMinimum(IntVec<3>& v, const IntVec<3> d, int* waveBz);

   template<>
   void shiftToMinimum(IntVec<2>& v, const IntVec<2> d, int* waveBz);

   template<>
   void shiftToMinimum(IntVec<1>& v, const IntVec<1> d, int* waveBz);
 

}
}
