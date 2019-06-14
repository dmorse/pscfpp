#ifndef PSCF_SPACE_GROUP_H
#define PSCF_SPACE_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/SpaceSymmetry.h>
#include <pscf/crystal/SymmetryGroup.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FSArray.h>
#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Group of crystal symmetries with no translations.
   *
   * \ingroup Crystal_Module
   */
   template <D>
   class SpaceGroup : public SymmetryGroup<SpaceSymmetry <D> >
   {

   public:

      /**
      * Generate a set of reciprocal vectors that are related by symmetry.
      *
      * \param root the reciprocal vector to which all others are related.
      * \param star array containing all members of the star (on ouput).
      */
      template <int N>
      void makeStar(const IntVec<D>& root, FSArray<IntVec<D>, N>& star);

   // friends:

      friend 
      std::ostream& operator << (std::ostream& out, const SpaceGroup& A);

   };

   /**
   * Output stream inserter operator writing a SpaceGroup to stream
   */ 
   std::ostream& operator << (std::ostream& out, const SpaceGroup& A);

   // Template function definition

   template <int N>
   void SpaceGroup::makeStar(const IntVec<D>& root, FSArray<IntVec<D>, N>& star)
   {
      // Precondition
      if (star.capacity() < size()) {
         UTIL_THROW("Star array capacity is less than order of group");
      }

      IntVec<D> vector;
      int       i, j;
      bool      found;
      star.clear();
      for (i = 0; i < size(); ++i) {
         vector = root*(*this)[i];
         found = false;
         for (j = 0; j < star.size(); ++j) {
            if (vector == star[j]) {
               found = true;
               break;
            }
         }
         if (!found) {
            star.append(vector);
         }
      }
   
   }

}
#endif
