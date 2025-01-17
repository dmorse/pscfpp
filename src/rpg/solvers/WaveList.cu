/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WaveList.tpp"

namespace Pscf { 
namespace Rpg 
{ 

   using namespace Util;

   // Explicit specializations
   template <>
   bool WaveList<1>::hasVariableAngle() const
   {  return false; }

   template <>
   bool WaveList<2>::hasVariableAngle() const
   {
      UTIL_CHECK(unitCell().lattice() != UnitCell<2>::Null);
      if ((unitCell().lattice() == UnitCell<2>::Oblique) ||
          (unitCell().lattice() == UnitCell<2>::Rhombic)) {
         return true;
      } else {
         return false;
      }
   }

   template <>
   bool WaveList<3>::hasVariableAngle() const
   {
      UTIL_CHECK(unitCell().lattice() != UnitCell<3>::Null);
      if ((unitCell().lattice() == UnitCell<3>::Monoclinic) ||
          (unitCell().lattice() == UnitCell<3>::Triclinic)  ||
          (unitCell().lattice() == UnitCell<3>::Rhombohedral)) {
         return true;
      } else {
         return false;
      }
   }

   template class WaveList<1>;
   template class WaveList<2>;
   template class WaveList<3>;

}
}
