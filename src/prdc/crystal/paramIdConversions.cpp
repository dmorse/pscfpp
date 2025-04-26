/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "paramIdConversions.h"

namespace Pscf {
namespace Prdc {

   template <> 
   int convertFullParamIdToReduced<1>(const int fullId, 
                  const typename UnitCell<1>::LatticeSystem lattice)
   {
      UTIL_CHECK(fullId == 0);
      UTIL_CHECK(lattice != UnitCell<1>::Null);
      return 0;
   }

   template <> 
   int convertFullParamIdToReduced<2>(const int fullId, 
                  const typename UnitCell<2>::LatticeSystem lattice)
   {
      UTIL_CHECK((fullId > -1) && (fullId < 3));
      UTIL_CHECK(lattice != UnitCell<2>::Null);

      if (fullId == 0) {
         return 0;
      } else if (fullId == 1) {
         if ((lattice == UnitCell<2>::Rectangular) ||
             (lattice == UnitCell<2>::Oblique)) {
            return 1;
         } else {
            return 0;
         } 
      } else { // fullId == 2
         if (lattice == UnitCell<2>::Oblique) {
            return 2;
         } else if (lattice == UnitCell<2>::Rhombic) {
            return 1;
         } else {
            return -1;
         } 
      }
   }

   template <> 
   int convertFullParamIdToReduced<3>(const int fullId, 
                  const typename UnitCell<3>::LatticeSystem lattice)
   {
      UTIL_CHECK((fullId > -1) && (fullId < 6));
      UTIL_CHECK(lattice != UnitCell<3>::Null);
      if (fullId == 0) {
         return 0;
      } else if (fullId == 1) { 
         if ((lattice == UnitCell<3>::Orthorhombic) ||
             (lattice == UnitCell<3>::Monoclinic) ||
             (lattice == UnitCell<3>::Triclinic)) {
            return 1;
         } else {
            return 0;
         }
      } else if (fullId == 2) { 
         if ((lattice == UnitCell<3>::Cubic) ||
             (lattice == UnitCell<3>::Rhombohedral)) {
            return 0;
         } else if ((lattice == UnitCell<3>::Tetragonal) ||
                    (lattice == UnitCell<3>::Hexagonal)) {
            return 1;
         } else { // lattice is Orthorhombic, Monoclinic, or Triclinic
            return 2;
         }
      } else if ((fullId == 3) || (fullId == 5)) {
         if (lattice == UnitCell<3>::Triclinic) {
            return fullId;
         } else if (lattice == UnitCell<3>::Rhombohedral) {
            return 1;
         } else {
            return -1;
         } 
      } else { // fullId == 4
         if (lattice == UnitCell<3>::Triclinic) {
            return fullId;
         } else if (lattice == UnitCell<3>::Monoclinic) {
            return 3;
         }else if (lattice == UnitCell<3>::Rhombohedral) {
            return 1;
         } else {
            return -1;
         } 
      }
   }

   template <>
   int convertReducedParamIdToFull<1>(const int reducedId, 
                  const typename UnitCell<1>::LatticeSystem lattice)
   {
      UTIL_CHECK(reducedId == 0);
      UTIL_CHECK(lattice != UnitCell<1>::Null);
      return 0;
   }

   template <>
   int convertReducedParamIdToFull<2>(const int reducedId, 
                  const typename UnitCell<2>::LatticeSystem lattice)
   {
      UTIL_CHECK(reducedId > -1);
      if ((lattice == UnitCell<2>::Square) ||
          (lattice == UnitCell<2>::Hexagonal)) {
         UTIL_CHECK(reducedId == 0);
         return reducedId;
      } else if (lattice == UnitCell<2>::Rectangular) {
         UTIL_CHECK(reducedId < 2);
         return reducedId;
      } else if (lattice == UnitCell<2>::Rhombic) {
         UTIL_CHECK(reducedId < 2);
         if (reducedId == 0) {
            return 0;
         } else { // reducedId == 1
            return 2;
         }
      } else { // lattice == Oblique
         UTIL_CHECK(reducedId < 3);
         return reducedId;
      }
   }

   template <>
   int convertReducedParamIdToFull<3>(const int reducedId, 
                  const typename UnitCell<3>::LatticeSystem lattice)
   {
      UTIL_CHECK(reducedId > -1);
      if (lattice == UnitCell<3>::Cubic) {
         UTIL_CHECK(reducedId == 0);
         return reducedId;
      } else if ((lattice == UnitCell<3>::Hexagonal) ||
                 (lattice == UnitCell<3>::Tetragonal)) {
         UTIL_CHECK(reducedId < 2);
         if (reducedId == 0) {
            return 0;
         } else { // reducedId == 1
            return 2;
         }
      } else if (lattice == UnitCell<3>::Rhombohedral) {
         UTIL_CHECK(reducedId < 2);
         if (reducedId == 0) {
            return 0;
         } else { // reducedId == 1
            return 3;
         }
      } else if (lattice == UnitCell<3>::Orthorhombic) {
         UTIL_CHECK(reducedId < 3);
         return reducedId;
      } else if (lattice == UnitCell<3>::Monoclinic) {
         UTIL_CHECK(reducedId < 4);
         if (reducedId < 3) {
            return reducedId;
         } else { // reducedId == 3
            return 4;
         }
      } else { // lattice == Triclinic
         UTIL_CHECK(reducedId < 6);
         return reducedId;
      }
   }

}
}
