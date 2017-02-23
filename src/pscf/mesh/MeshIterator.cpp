/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MeshIterator.h"

namespace Pscf
{

   template <>
   void MeshIterator<1>::operator ++ ()
   {
      position_[0]++;
      if (position_[0] == dimensions_[0]) {
         position_[0] = 0;
      }
      rank_++;
   }

   template <>
   inline void MeshIterator<2>::operator ++ ()
   {
      position_[1]++;
      if (position_[1] == dimensions_[1]) {
         position_[1] = 0;
         increment(0);
      }
      rank_++;
   }

   template <>
   inline void MeshIterator<3>::operator ++ ()
   {
      position_[2]++;
      if (position_[2] == dimensions_[2]) {
         position_[2] = 0;
         increment(1);
      }
      rank_++;
   }

}
