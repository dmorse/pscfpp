#ifndef PSCF_SPACE_GROUP_TPP
#define PSCF_SPACE_GROUP_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/SpaceGroup.h>
#include <pscf/crystal/groupFile.h>

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
   const
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

   template <int D>
   void SpaceGroup<D>::shiftOrigin(
                    typename SpaceSymmetry<D>::Translation const & origin)
   {
      for (int i = 0; i < size(); ++i) {
         (*this)[i].shiftOrigin(origin);
      }
   }

   /*
   * Read a group from file
   */
   template <int D>
   void readGroup(std::string groupName, SpaceGroup<D>& group)
   {
      if (groupName == "I") {
         // Create identity group by default
         group.makeCompleteGroup();
      } else {
         bool foundFile = false;
         {
            // Search first in this directory
            std::ifstream in;
            in.open(groupName);
            if (in.is_open()) {
               in >> group;
               UTIL_CHECK(group.isValid());
               foundFile = true;
            }
         }
         if (!foundFile) {
            // Search in the data directory containing standard space groups
            std::string fileName = makeGroupFileName(D, groupName);
            std::ifstream in;
            in.open(fileName);
            if (in.is_open()) {
               in >> group;
               UTIL_CHECK(group.isValid());
            } else {
               Log::file() << "\nFailed to open group file: " 
                           << fileName << "\n";
               Log::file() << "\n Error: Unknown space group\n";
               UTIL_THROW("Unknown space group");
            }
         } 
      }
   }

   /*
   * Open an output file and write group to file. 
   */
   template <int D>
   void writeGroup(std::string filename, SpaceGroup<D> const & group)
   {
      std::ofstream out;
      out.open(filename);
      out << group;
   }

}
#endif
