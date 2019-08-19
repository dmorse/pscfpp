/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
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
   UnitCell<1>::UnitCell()
    : lattice_(Null)
   {}

   /*
   * Set nParameter based on the lattice system.
   */
   void UnitCell<1>::setNParameter()
   {
      UTIL_CHECK(lattice_ != UnitCell<1>::Null);
      if (lattice_ == UnitCell<1>::Lamellar) {
         nParameter_ = 1;
      } else {
         UTIL_THROW("Invalid lattice system value");
      }
   }

   /*
   * Set the Bravais and reciprocal lattice parameters.
   */
   void UnitCell<1>::setBasis()
   {
      UTIL_CHECK(lattice_ != UnitCell<1>::Null);
      UTIL_CHECK(nParameter_ == 1);

      rBasis_[0][0] = parameters_[0];
      kBasis_[0][0] = 2.0*Constants::Pi/parameters_[0];
      drBasis_[0](0,0) = 1.0;
   }

   /*
   * Extract a UnitCell<1>::LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<1>::LatticeSystem& lattice)
   {

      std::string buffer;
      in >> buffer;
      if (buffer == "Lamellar" || buffer == "lamellar") {
         lattice = UnitCell<1>::Lamellar;
      } else {
         UTIL_THROW("Invalid UnitCell<1>::LatticeSystem value input");
      }
      return in;
   }

   /*
   * Insert a UnitCell<1>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<1>::LatticeSystem lattice)
   {
      if (lattice == UnitCell<1>::Lamellar) {
         out << "lamellar";
      } else {
         UTIL_THROW("Invalid value of UnitCell<1>::Lamellar");
      }
      return out;
   }

}
