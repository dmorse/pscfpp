/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"
#include <util/math/Constants.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor.
   */
   UnitCell<2>::UnitCell()
    : lattice_(Null)
   {}

   /*
   * Read the lattice system and set nParameter.
   */
   void UnitCell<2>::setNParameter()
   {
      UTIL_CHECK(lattice_ != UnitCell<2>::Null);
      if (lattice_ == UnitCell<2>::Square) {
         nParameter_ = 1;
      } else
      if (lattice_ == UnitCell<2>::Hexagonal) {
         nParameter_ = 1;
      } else
      if (lattice_ == UnitCell<2>::Rectangular) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<2>::Rhombic) {
         nParameter_ = 2;
      } else
      if (lattice_ == UnitCell<2>::Oblique) {
         nParameter_ = 3;
      } else {
         UTIL_THROW("Invalid lattice system value");
      }
   }

   /*
   * Set the Bravais and reciprocal lattice vectors.
   */
   void UnitCell<2>::setBasis()
   {
      UTIL_CHECK(lattice_ != UnitCell<2>::Null);
      UTIL_CHECK(nParameter_ > 0);

      double twoPi = 2.0*Constants::Pi;
      int i;
      if (lattice_ == UnitCell<2>::Square) {
         UTIL_CHECK(nParameter_ == 1);
         for (i=0; i < 2; ++i) { 
            rBasis_[i][i] = parameters_[0];
            kBasis_[i][i] = twoPi/parameters_[0];
            drBasis_[0](i,i) = 1.0;
         }
      } else 
      if (lattice_ == UnitCell<2>::Rectangular) {
         UTIL_CHECK(nParameter_ == 2);
         for (i=0; i < 2; ++i) { 
            rBasis_[i][i] = parameters_[i];
            kBasis_[i][i] = twoPi/parameters_[i];
            drBasis_[i](i,i) = 1.0;
         }
      } else 
      if (lattice_ == UnitCell<2>::Hexagonal) {
         UTIL_CHECK(nParameter_ == 1);
         double a = parameters_[0];
         double rt3 = sqrt(3.0);

         rBasis_[0][0] = a;
         rBasis_[0][1] = 0.0;
         rBasis_[1][0] = -0.5*a;
         rBasis_[1][1] = 0.5*rt3*a;

         drBasis_[0](0, 0) = 1.0;
         drBasis_[0](0, 1) = 0.0;
         drBasis_[0](1, 0) = -0.5;
         drBasis_[0](1, 1) = 0.5*rt3;

         kBasis_[0][0] = twoPi/a;
         kBasis_[0][1] = twoPi/(rt3*a);
         kBasis_[1][0] = 0.0;
         kBasis_[1][1] = twoPi/(0.5*rt3*a);
      } else 
      if (lattice_ == UnitCell<2>::Rhombic) {
         UTIL_CHECK(nParameter_ == 2);
         double a = parameters_[0];
         double gamma = parameters_[1];
         // gamma is the angle between the two Bravais basis vectors

         double cg = cos(gamma);
         double sg = sin(gamma);

         
         rBasis_[0][0] = a;
         rBasis_[0][1] = 0.0;
         rBasis_[1][0] = cg*a;
         rBasis_[1][1] = sg*a;

         drBasis_[0](0, 0) = 1.0;
         drBasis_[0](0, 1) = 0.0;
         drBasis_[0](1, 0) = cg;
         drBasis_[0](1, 1) = sg;
         drBasis_[1](1, 0) = -sg*a;
         drBasis_[1](1, 1) =  cg*a;

         kBasis_[0][0] = twoPi/a;
         kBasis_[0][1] = -twoPi*cg/(sg*a);
         kBasis_[1][0] = 0.0;
         kBasis_[1][1] = twoPi/(a*sg);

      } else 
      if (lattice_ == UnitCell<2>::Oblique) {
         UTIL_CHECK(nParameter_ == 3);
         double a = parameters_[0];
         double b = parameters_[1];
         double gamma = parameters_[2];
         // gamma is the angle between the two Bravais basis vectors

         double cg = cos(gamma);
         double sg = sin(gamma);
         
         rBasis_[0][0] = a;
         rBasis_[0][1] = 0.0;
         rBasis_[1][0] = cg*b;
         rBasis_[1][1] = sg*b;

         drBasis_[0](0, 0) = 1.0;
         drBasis_[0](0, 1) = 0.0;
         drBasis_[1](1, 0) = cg;
         drBasis_[1](1, 1) = sg;
         drBasis_[1](1, 0) = -sg*b;
         drBasis_[1](1, 1) =  cg*b;

         kBasis_[0][0] = twoPi/a;
         kBasis_[0][1] = -twoPi*cg/(sg*a);
         kBasis_[1][0] = 0.0;
         kBasis_[1][1] = twoPi/(b*sg);

      } else {
         UTIL_THROW("Unimplemented 2D lattice type");
      }
   }

   /*
   * Extract a UnitCell<2>::LatticeSystem from an istream as a string.
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<2>::LatticeSystem& lattice)
   {

      std::string buffer;
      in >> buffer;
      if (buffer == "Square" || buffer == "square") {
         lattice = UnitCell<2>::Square;
      } else
      if (buffer == "Rectangular" || buffer == "rectangular") {
         lattice = UnitCell<2>::Rectangular;
      } else
      if (buffer == "Rhombic" || buffer == "rhombic") {
         lattice = UnitCell<2>::Rhombic;
      } else
      if (buffer == "Hexagonal" || buffer == "hexagonal") {
         lattice = UnitCell<2>::Hexagonal;
      } else
      if (buffer == "Oblique" || buffer == "oblique") {
         lattice = UnitCell<2>::Oblique;
      } else {
         UTIL_THROW("Invalid UnitCell<2>::LatticeSystem value input");
      }
      return in;
   }

   /*
   * Insert a UnitCell<2>::LatticeSystem to an ostream as a string.
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<2>::LatticeSystem lattice)
   {
      if (lattice == UnitCell<2>::Square) {
         out << "square";
      } else
      if (lattice == UnitCell<2>::Rectangular) {
         out << "rectangular";
      } else
      if (lattice == UnitCell<2>::Rhombic) {
         out << "rhombic";
      } else
      if (lattice == UnitCell<2>::Hexagonal) {
         out << "hexagonal";
      } else
      if (lattice == UnitCell<2>::Oblique) {
         out << "oblique";
      } else
      if (lattice == UnitCell<2>::Null) {
         out << "Null";
      } else {
         UTIL_THROW("This should never happen");
      }
      return out;
   }

   /*
   * Assignment operator.
   */
   UnitCell<2>& UnitCell<2>::operator = (const UnitCell<2>& other)
   {
      if (lattice_ != UnitCell<2>::Null) {
         UTIL_CHECK(other.lattice_ == lattice_);
      }
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
   * Set state of the unit cell (lattice system and parameters)
   */
   void UnitCell<2>::set(UnitCell<2>::LatticeSystem lattice)
   {
      UTIL_CHECK(lattice != UnitCell<2>::Null);
      if (lattice_ != UnitCell<2>::Null) {
         UTIL_CHECK(lattice == lattice_);
      }
      isInitialized_ = false;
      lattice_ = lattice;
      setNParameter();
   }

   /*
   * Set state of the unit cell (lattice system and parameters)
   */
   void UnitCell<2>::set(UnitCell<2>::LatticeSystem lattice,
                         FSArray<double, 6> const & parameters)
   {
      set(lattice);
      UTIL_CHECK(parameters.size() == nParameter_);
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i] = parameters[i];
      }   
      setLattice();
   }

   /*
   * Get the generalized 2D volume (area) of the unit cell. 
   */
   double UnitCell<2>::volume() const
   {
      double a = 0.0;
      a += rBasis_[0][0]*rBasis_[1][1];
      a -= rBasis_[0][1]*rBasis_[1][0];
      return a;
   }

}
}
