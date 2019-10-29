/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"
#include <util/math/Constants.h>

namespace Pscf
{

   using namespace Util;

   /*
   * Constructor.
   */
   UnitCell<3>::UnitCell()
    : lattice_(Null)
   {}

   /*
   * Read the lattice system and set nParameter.
   */
   void UnitCell<3>::setNParameter()
   {
      UTIL_CHECK(lattice_ != UnitCell<3>::Null);
      if (lattice_ == UnitCell<3>::Cubic) {
         nParameter_ = 1;
      } else
      if (lattice_ == UnitCell<3>::Tetragonal) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<3>::Orthorhombic) {
         nParameter_ = 3;
      } else
      if (lattice_ == UnitCell<3>::Monoclinic) {
         nParameter_ = 4;
      } else
      if (lattice_ == UnitCell<3>::Triclinic) {
         nParameter_ = 6;
      } else
      if (lattice_ == UnitCell<3>::Rhombohedral) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<3>::Hexagonal) {
         nParameter_ = 2;
      } else {
         UTIL_THROW("Invalid value");
      }
   }

   /*
   * Set the Bravais and reciprocal lattice parameters.
   */
   void UnitCell<3>::setBasis()
   {
      UTIL_CHECK(lattice_ != UnitCell<3>::Null);
      UTIL_CHECK(nParameter_ > 0);

      // Set elements for specific lattice types
      double twoPi = 2.0*Constants::Pi;
      int i;
      if (lattice_ == UnitCell<3>::Cubic) {
         UTIL_CHECK(nParameter_ == 1);
         for (i=0; i < 3; ++i) { 
            rBasis_[i][i] = parameters_[0];
            kBasis_[i][i] = twoPi/parameters_[0];
            drBasis_[0](i,i) = 1.0;
         }
      } else 
      if (lattice_ == UnitCell<3>::Tetragonal) {
         UTIL_CHECK(nParameter_ == 2);
         rBasis_[0][0] = parameters_[0];
         rBasis_[1][1] = rBasis_[0][0];
         rBasis_[2][2] = parameters_[1];
         kBasis_[0][0] = twoPi/parameters_[0];
         kBasis_[1][1] = kBasis_[0][0];
         kBasis_[2][2] = twoPi/parameters_[1];
         drBasis_[0](0,0) = 1.0;
         drBasis_[0](1,1) = 1.0;
         drBasis_[1](2,2) = 1.0;
      } else 
      if (lattice_ == UnitCell<3>::Hexagonal) {
         UTIL_CHECK(nParameter_ == 2);
         double a = parameters_[0];
         double c = parameters_[1];
         double rt3 = sqrt(3.0);

         rBasis_[0][0] = a;
         rBasis_[0][1] = 0.0;
         rBasis_[1][0] = -0.5*a;
         rBasis_[1][1] = 0.5*rt3*a;
         rBasis_[2][2] = c;

         drBasis_[0](0, 0) = 1.0;
         drBasis_[0](0, 1) = 0.0;
         drBasis_[0](1, 0) = -0.5;
         drBasis_[0](1, 1) = 0.5*rt3;
         drBasis_[1](2, 2) = 1.0;

         kBasis_[0][0] = twoPi/a;
         kBasis_[0][1] = twoPi/(rt3*a);
         kBasis_[1][0] = 0.0;
         kBasis_[1][1] = twoPi/(0.5*rt3*a);
         kBasis_[2][2] = twoPi/(c);
 
      }else 
      if (lattice_ == UnitCell<3>::Orthorhombic) {
         UTIL_CHECK(nParameter_ == 3);
         for (i=0; i < 3; ++i) { 
            rBasis_[i][i] = parameters_[i];
            kBasis_[i][i] = twoPi/parameters_[i];
            drBasis_[i](i,i) = 1.0;
         }
      } else {
         UTIL_THROW("Unimplemented 3D lattice type");
      }
   }

   /*
   * Extract a UnitCell<3>::LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<3>::LatticeSystem& lattice)
   {

      std::string buffer;
      in >> buffer;
      if (buffer == "Cubic" || buffer == "cubic") {
         lattice = UnitCell<3>::Cubic;
      } else
      if (buffer == "Tetragonal" || buffer == "tetragonal") {
         lattice = UnitCell<3>::Tetragonal;
      } else
      if (buffer == "Orthorhombic" || buffer == "orthorhombic") {
         lattice = UnitCell<3>::Orthorhombic;
      } else
      if (buffer == "Monoclinic" || buffer == "monoclinic") {
         lattice = UnitCell<3>::Monoclinic;
      } else
      if (buffer == "Triclinic" || buffer == "triclinic") {
         lattice = UnitCell<3>::Triclinic;
      } else
      if (buffer == "Rhombohedral" || buffer == "rhombohedral") {
         lattice = UnitCell<3>::Rhombohedral;
      } else
      if (buffer == "Hexagonal" || buffer == "hexagonal") {
         lattice = UnitCell<3>::Hexagonal;
      } else {
         UTIL_THROW("Invalid UnitCell<3>::LatticeSystem value input");
      }
      return in;
   }

   /*
   * Insert a UnitCell<3>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<3>::LatticeSystem lattice)
   {
      if (lattice == UnitCell<3>::Cubic) {
         out << "cubic";
      } else
      if (lattice == UnitCell<3>::Tetragonal) {
         out << "tetragonal";
      } else
      if (lattice == UnitCell<3>::Orthorhombic) {
         out << "orthorhombic";
      } else
      if (lattice == UnitCell<3>::Monoclinic) {
         out << "monoclinic";
      } else
      if (lattice == UnitCell<3>::Triclinic) {
         out << "triclinic";
      } else
      if (lattice == UnitCell<3>::Rhombohedral) {
         out << "rhombohedral";
      } else
      if (lattice == UnitCell<3>::Hexagonal) {
         out << "hexagonal";
      } else {
         UTIL_THROW("This should never happen");
      }
      return out;
   }

}
