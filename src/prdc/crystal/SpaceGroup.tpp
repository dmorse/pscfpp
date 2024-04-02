#ifndef PRDC_SPACE_GROUP_TPP
#define PRDC_SPACE_GROUP_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/groupFile.h>

namespace Pscf {
namespace Prdc {

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
   * Check if a mesh is compatible with this space group.
   */
   template <int D>
   void SpaceGroup<D>::checkMeshDimensions(IntVec<D> const & dimensions)
   const
   {

      // ---------------------------------------------------------------
      // Check compatibility of mesh dimensions & translations
      // ---------------------------------------------------------------

      // Identify required divisor of mesh dimension in each direction
      int numerator, denominator;
      IntVec<D> divisors;
      for (int i = 0; i < D; ++i) {
         // Find maximum denominator
         divisors[i] = 1;
         for (int j = 0; j < size(); ++j) {
            numerator = (*this)[j].t(i).num();
            denominator = (*this)[j].t(i).den();
            if (numerator != 0) {
               UTIL_CHECK(denominator > 0);
               if (denominator > divisors[i]) {
                  divisors[i] = denominator;
               } 
            }
         }
         // Make sure each divisor is a multiple of all denominators
         for (int j = 0; j < size(); ++j) {
            numerator = (*this)[j].t(i).num();
            denominator = (*this)[j].t(i).den();
            if (numerator != 0) {
               if (denominator < divisors[i]) {
                  if (divisors[i]%denominator != 0) {
                     divisors[i] = divisors[i]*denominator;
                  }
               }
            }
         }
      }

      // Check that mesh dimensions are multiples of required divisor
      for (int i = 1; i < D; ++i) {
         if (dimensions[i]%divisors[i] != 0) {
            Log::file() 
               << "\n"
               << "Mesh dimensions incompatible with the space group:\n" 
               << "  dimensions[" << i << "] = " << dimensions[i] << "\n"
               << "  This dimension must be a multiple of " << divisors[i] 
               << "\n\n";
               UTIL_THROW("Error: Mesh not incompatible with space group");
         }
      }

      // ---------------------------------------------------------------
      // Check compatibility of mesh dimensions & point group operations
      // ---------------------------------------------------------------

      // Identify pairs of directions that are related by point group ops
      FMatrix<bool, D, D> areRelated;
      for (int i = 0; i < D; ++i) {
         for (int j = 0; j < D; ++j) {
            areRelated(i, j) = false;
         }
         areRelated(i, i) = true;
      }
      for (int k = 0; k < size(); ++k) {
         for (int i = 0; i < D; ++i) {
            for (int j = 0; j < D; ++j) {
               if (i != j) {
                  if ( (*this)[k].R(i,j) != 0) {
                     areRelated(i, j) = true;
                     areRelated(j, i) = true;
                  }
               }
            }
         }
      }

      // Check if mesh dimensions are equal for related directions
      for (int i = 0; i < D; ++i) {
         for (int j = 0; j < i; ++j) {
            if (areRelated(i,j) && (dimensions[i] != dimensions[j])) {
               Log::file() 
                 << "\n"
                 << "Mesh dimensions incompatible with the space group - \n"
                 << "Unequal dimensions for related directions:\n"
                 << "  dimensions[" << i << "] = " << dimensions[i] << "\n"
                 << "  dimensions[" << j << "] = " << dimensions[j] 
                 << "\n\n";
               UTIL_THROW("Error: Mesh not incompatible with space group");
            }
         }
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
}
#endif
