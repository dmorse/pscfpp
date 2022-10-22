/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
      } else 
      if (lattice == UnitCell<1>::Null) {
         out << "Null";
      } else {
         UTIL_THROW("Invalid value of UnitCell<1>::Lamellar");
      }
      return out;
   }

   /*
   * Assignment operator.
   */
   UnitCell<1>& UnitCell<1>::operator = (const UnitCell<1>& other)
   {
      isInitialized_ = false;
      lattice_ = other.lattice_;
      setNParameter();
      UTIL_CHECK(nParameter_ == other.nParameter_);
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i] = other.parameters_[i];
      }   
      setLattice();
      return *this;
   }

   /*
   * Set the lattice system, but not unit cell parameters.
   */
   void UnitCell<1>::set(UnitCell<1>::LatticeSystem lattice)
   {
      isInitialized_ = false;
      lattice_ = lattice;
      setNParameter();
   }

   /*
   * Set state of the unit cell (lattice system and parameters).
   */
   void UnitCell<1>::set(UnitCell<1>::LatticeSystem lattice,
                         FSArray<double, 6> const & parameters)
   {
      isInitialized_ = false;
      lattice_ = lattice;
      setNParameter();
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i] = parameters[i];
      }   
      setLattice();
   }

}
