#ifndef UTIL_POINT_GROUP_H
#define UTIL_POINT_GROUP_H

#include <util/crystal/PointSymmetry.h>
#include <util/crystal/SymmetryGroup.h>
#include <util/containers/FSArray.h>
#include <util/space/IntVector.h>

#include <iostream>

namespace Util {

   /**
   * Group of crystal symmetries with no translations.
   *
   * \ingroup Crystal_Module
   */
   class PointGroup : public SymmetryGroup<PointSymmetry>
   {

   public:

      /**
      * Generate a set of reciprocal vectors that are related by symmetry.
      *
      * \param root the reciprocal vector to which all others are related.
      * \param star array containing all members of the star (on ouput).
      */
      template <int N>
      void makeStar(const IntVector& root, FSArray<IntVector, N>& star);

   // friends:

      friend 
      std::ostream& operator << (std::ostream& out, const PointGroup& A);

   };

   /**
   * Output stream inserter operator writing a PointGroup to stream
   */ 
   std::ostream& operator << (std::ostream& out, const PointGroup& A);

   // Template function definition

   template <int N>
   void PointGroup::makeStar(const IntVector& root, FSArray<IntVector, N>& star)
   {
      // Precondition
      if (star.capacity() < size()) {
         UTIL_THROW("Star array capacity is less than order of group");
      }

      IntVector vector;
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
