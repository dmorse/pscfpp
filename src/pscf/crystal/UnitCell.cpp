/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"

namespace Pscf
{ 

   using namespace Util;

   template class UnitCell<1>;
   template class UnitCell<2>;
   template class UnitCell<3>;

   #if 0
   /* 
   * Extract a UnitCell<3>::LatticeSystem from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, 
                            UnitCell<3>::LatticeSystem& lattice)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Cubic" || buffer == "cubic") {
         lattice = Cubic;
      } else 
      if (buffer == "Tetragonal" || buffer == "tetragonal") {
         lattice = Tetragonal;
      } else 
      if (buffer == "Orthorhombic" || buffer == "orthorhombic") {
         lattice = Orthorhombic;
      } else 
      if (buffer == "Monoclinic" || buffer == "monoclinic") {
         lattice = Monoclinic; 
      } else 
      if (buffer == "Triclinic" || buffer == "triclinic") {
         lattice = Triclinic; 
      } else 
      if (buffer == "Rhombohedral" || buffer == "rhombohedral") {
         lattice = Rhombohedral;
      } else 
      if (buffer == "Hexagonal" || buffer == "hexagonal") {
         lattice = Hexagonal; 
      } else {
         #ifndef UTIL_MPI
         Log::file() << "Unknown UnitCell<3>::LatticeSystem: " 
                     << buffer << std::endl;
         #endif
         UTIL_THROW("Invalid UnitCell<3>::LatticeSystem value input");
      }
      return in;
   }
   
   /* 
   * Insert a UnitCell<3>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, 
                            UnitCell<3>::LatticeSystem lattice) 
   {
      if (lattice == Cubic) {
         out << "cubic";
      } else 
      if (lattice == Tetragonal) {
         out << "tetragonal";
      } else
      if (lattice == Orthorhombic) {
         out << "orthorhombic";
      } else
      if (lattice == Monoclinic) {
         out << "monoclinic";
      } else
      if (lattice == Triclinic) {
         out << "triclinic";
      } else
      if (lattice == Rhombohedral) {
         out << "rhombohedral";
      } else
      if (lattice == Hexagonal) {
         out << "hexagonal";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out; 
   }
   #endif

} 
