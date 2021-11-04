#ifndef PSCF_SPACE_GROUP_TPP
#define PSCF_SPACE_GROUP_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/SpaceGroup.h>

namespace Pscf
{

   using namespace Util;

   /*
   * Find inversion center, if any.
   */
   template <int D>
   bool 
   SpaceGroup<D>::hasInversionCenter(
                            typename SpaceSymmetry<D>::Translation& center)
   {
      bool isInversion = false;
      int i, j, k;
      for (i = 0; i < size(); ++i) {
         isInversion = true;
         for (j = 0; j < D; ++j) {
            for (k = 0; k < D; ++k) {
               if (j == k) {
                  if ((*this)[i].R(j,k) != -1) isInversion = false;
               } else {
                  if ((*this)[i].R(j,k) !=  0) isInversion = false;
               }
            }
         }
         if (isInversion) {
            for (int j = 0; j < D; ++j) {
               center[j] = (*this)[i].t(j)/2;
            }
            return true;
         }
      }
      return false;
   }

}
#endif
