/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/hasVariableAngle.h>

namespace Pscf { 
namespace Prdc { 

   template <>
   bool hasVariableAngle<1>(typename UnitCell<1>::LatticeSystem lattice)
   {  return false; }

   template <>
   bool hasVariableAngle<2>(typename UnitCell<2>::LatticeSystem lattice)
   {
      UTIL_CHECK(lattice != UnitCell<2>::Null);
      if ((lattice == UnitCell<2>::Oblique) ||
          (lattice == UnitCell<2>::Rhombic)) {
         return true;
      } else {
         return false;
      }
   }

   template <>
   bool hasVariableAngle<3>(typename UnitCell<3>::LatticeSystem lattice)
   {
      UTIL_CHECK(lattice != UnitCell<3>::Null);
      if ((lattice == UnitCell<3>::Monoclinic) ||
          (lattice == UnitCell<3>::Triclinic)  ||
          (lattice == UnitCell<3>::Rhombohedral)) {
         return true;
      } else {
         return false;
      }
   }

}
}
