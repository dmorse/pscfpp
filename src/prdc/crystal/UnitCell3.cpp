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
   * Set Bravais and reciprocal lattice parameters, and drBasis.
   * 
   * Function UnitCellBase::initializeToZero() must be called before
   * setBasis() to initialize all elements of vectors and tensors to
   * zero. Only nonzero elements need to be set here.
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
         rBasis_[1][1] = parameters_[0];
         rBasis_[2][2] = parameters_[1];

         drBasis_[0](0,0) = 1.0;
         drBasis_[0](1,1) = 1.0;
         drBasis_[1](2,2) = 1.0;

         kBasis_[0][0] = twoPi/parameters_[0];
         kBasis_[1][1] = kBasis_[0][0];
         kBasis_[2][2] = twoPi/parameters_[1];

      }else 
      if (lattice_ == UnitCell<3>::Orthorhombic) {
         UTIL_CHECK(nParameter_ == 3);
         for (i=0; i < 3; ++i) { 
            rBasis_[i][i] = parameters_[i];
            drBasis_[i](i,i) = 1.0;
            kBasis_[i][i] = twoPi/parameters_[i];
         }
      } else 
      if (lattice_ == UnitCell<3>::Hexagonal) {
         UTIL_CHECK(nParameter_ == 2);
         double a = parameters_[0];
         double c = parameters_[1];
         double rt3 = sqrt(3.0);

         rBasis_[0][0] = a;
         rBasis_[1][0] = -0.5*a;
         rBasis_[1][1] = 0.5*rt3*a;
         rBasis_[2][2] = c;

         drBasis_[0](0, 0) = 1.0;
         drBasis_[0](1, 0) = -0.5;
         drBasis_[0](1, 1) = 0.5*rt3;
         drBasis_[1](2, 2) = 1.0;

         kBasis_[0][0] = twoPi/a;
         kBasis_[0][1] = twoPi/(rt3*a);
         kBasis_[1][0] = 0.0;
         kBasis_[1][1] = twoPi/(0.5*rt3*a);
         kBasis_[2][2] = twoPi/(c);

      } else 
      if (lattice_ == UnitCell<3>::Rhombohedral) {
         UTIL_CHECK(nParameter_ == 2);

         // Set parameters
         double a, beta, theta;
         a = parameters_[0];    // magnitude of Bravais basis vectors
         beta = parameters_[1]; // angle between basis vectors
         theta = acos( sqrt( (2.0*cos(beta) + 1.0)/3.0 ) );
         // theta is the angle of all basis vectors from the z axis

         /*
         * Bravais lattice vectora a_i with i=0, 1, or 2 have form:
         *
         *    a_i = a * (z * cos(theta) + u_i * sin(theta) )
         * 
         * where z denotes a unit vector parallel to the z axis, and 
         * the u_i are three unit vectors in the x-y separated by 
         * angles of 2 pi /3  (or 120 degrees), for which u_0 is
         * parallel to the x axis. Note that u_i . u_j = -1/2 for any
         * i not equal to j.
         * 
         * The angle between any two such Bravais basis vectors is easily
         * computed from the dot product, a_i . a_j for unequal i and j.
         * This gives:
         * 
         *    cos beta = cos^2(theta) - 0.5 sin^2 (theta)
         *             = (3cos^{2}(theta) - 1)/2
         *
         * The angle theta is computed above by solving this equation
         * for theta as a function of the input parameter beta.
         *
         * The corresponding reciprocal lattice vectors have the form
         *
         *   b_i = (2*pi/3*a)( z/cos(theta) + 2 * u_i / sin(theta) )
         *
         * It is straightforward to confirm that these reciprocal
         * vectors satisfy the condition a_i . b_j = 2 pi delta_{ij}
         */

         // Trig function values
         double cosT = cos(theta);
         double sinT = sin(theta);
         double sqrt3 = sqrt(3.0);
         double sinPi3 = 0.5*sqrt3;
         double cosPi3 = 0.5;

         // Bravais basis vectors
         rBasis_[0][0] = sinT * a;
         rBasis_[0][2] = cosT * a;
         rBasis_[1][0] = -cosPi3 * sinT * a;
         rBasis_[1][1] = sinPi3 * sinT * a;
         rBasis_[1][2] = cosT * a;
         rBasis_[2][0] = -cosPi3 * sinT * a;
         rBasis_[2][1] = -sinPi3 * sinT * a;
         rBasis_[2][2] = cosT * a;

         // Reciprocal lattice basis vectors 
         kBasis_[0][0] = 2.0 * twoPi /( 3.0 * sinT * a);
         kBasis_[0][2] = twoPi /( 3.0 * cosT * a);
         kBasis_[1][0] = -2.0 * cosPi3 * twoPi /( 3.0 * sinT * a);
         kBasis_[1][1] =  2.0 * sinPi3 * twoPi /( 3.0 * sinT * a);
         kBasis_[1][2] = twoPi /( 3.0 * cosT * a);
         kBasis_[2][0] = -2.0 * cosPi3 * twoPi /( 3.0 * sinT * a);
         kBasis_[2][1] = -2.0 * sinPi3 * twoPi /( 3.0 * sinT * a);
         kBasis_[2][2] = twoPi /( 3.0 * cosT * a);

         // Derivatives with respect to length a
         drBasis_[0](0,0) = sinT;
         drBasis_[0](0,2) = cosT;
         drBasis_[0](1,0) = -cosPi3 * sinT;
         drBasis_[0](1,1) = sinPi3 * sinT;
         drBasis_[0](1,2) = cosT;
         drBasis_[0](2,0) = -cosPi3 * sinT;
         drBasis_[0](2,1) = -sinPi3* sinT;
         drBasis_[0](2,2) = cosT;

         // Define alpha = d(theta)/d(beta)
         double alpha = -sin(beta)/(3.0*cosT*sinT);

         // Derivatives with respect to beta
         drBasis_[1](0,0) = cosT * alpha * a;
         drBasis_[1](0,2) = -sinT * alpha * a;
         drBasis_[1](1,0) = -cosPi3 * cosT * alpha * a;
         drBasis_[1](1,1) = sinPi3 * cosT * alpha * a;
         drBasis_[1](1,2) = -sinT * alpha * a;
         drBasis_[1](2,0) = -cosPi3 * cosT * alpha * a;
         drBasis_[1](2,1) = -sinPi3* cosT * alpha * a;
         drBasis_[1](2,2) = -sinT * alpha * a;
      
      } else 
      if (lattice_ == UnitCell<3>::Monoclinic) {
         UTIL_CHECK(nParameter_ == 4);

         /*
         * Description: Bravais basis vectors A, B, and C
         * A (vector 0) is along x axis
         * B (vector 1) is along y axis
         * C (vector 2) is in the x-z plane, an angle beta from A/x
         * B is the unique axis (perpendicular to A and C)
         * Angle beta is stored in radians.
         */

         // Parameters
         double a = parameters_[0];     // length of A
         double b = parameters_[1];     // length of B
         double c = parameters_[2];     // length of C
         double beta  = parameters_[3]; // angle between C and A/x

         double cb = cos(beta);
         double sb = sin(beta);

         // Nonzero components of basis vectors
         // For rBasis_[i][j], i:basis vector, j:Cartesian component
         rBasis_[0][0] = a;
         rBasis_[1][1] = b;
         rBasis_[2][0] = c*cb;
         rBasis_[2][2] = c*sb;

         // Nonzero derivatives of basis vectors with respect to parameters
         // For drBasis_[k](i,j), k:parameter, i:basis vector, j:Cartesian
         drBasis_[0](0,0) = 1.0;
         drBasis_[1](1,1) = 1.0;
         drBasis_[2](2,0) = cb;
         drBasis_[2](2,2) = sb;
         drBasis_[3](2,0) = c*sb;
         drBasis_[3](2,2) = -c*cb;

         // Reciprocal lattice vectors
         kBasis_[0][0] = twoPi/a;
         kBasis_[0][2] = -twoPi*cb/(a*sb);
         kBasis_[1][1] = twoPi / b;
         kBasis_[2][2] = twoPi/(c*sb);

      } else 
      if (lattice_ == UnitCell<3>::Triclinic) {
         UTIL_CHECK(nParameter_ == 6);

         /*
         * Description: Bravais basis vectors A, B, and C
         * A (vector 0) is along x axis
         * B (vector 1) is in the x-y plane, an angle gamma from x axis
         * C (vector 2) is tilted by an angle theta from z axis
         * phi is the angle beween c-z an x-z (or a-z) planes
         * All three angles are stored in radians
         */
         double a = parameters_[0];      // length of A
         double b = parameters_[1];      // length of B
         double c = parameters_[2];      // length of C
         double phi = parameters_[3];    // angle between c-z and a-z planes
         double theta  = parameters_[4]; // angle between c and z axis
         double gamma = parameters_[5];  // angle between a and b

         // sine and cosine of all angles
         double cosPhi = cos(phi);
         double sinPhi = sin(phi);
         double cosTheta = cos(theta);
         double sinTheta = sin(theta);
         double cosGamma = cos(gamma);
         double sinGamma = sin(gamma);

         // Nonzero components of Bravais basis vectors
         rBasis_[0][0] = a;
         rBasis_[1][0] = b * cosGamma;
         rBasis_[1][1] = b * sinGamma;
         rBasis_[2][0] = c * sinTheta * cosPhi;
         rBasis_[2][1] = c * sinTheta * sinPhi;
         rBasis_[2][2] = c * cosTheta;

         // Nonzero derivatives of basis vectors with respect to parameters
         drBasis_[0](0,0) = 1.0;
         drBasis_[1](1,0) = cosGamma;
         drBasis_[1](1,1) = sinGamma;
         drBasis_[2](2,0) = sinTheta * cosPhi;
         drBasis_[2](2,1) = sinTheta * sinPhi;
         drBasis_[2](2,2) = cosTheta;
         drBasis_[3](2,0) = -c * sinTheta * sinPhi;
         drBasis_[3](2,1) =  c * sinTheta * cosPhi;
         drBasis_[4](2,0) =  c * cosTheta * cosPhi;
         drBasis_[4](2,1) =  c * cosTheta * sinPhi;
         drBasis_[4](2,2) = -c * sinTheta;
         drBasis_[5](1,0) = -b * sinGamma;
         drBasis_[5](1,1) =  b * cosGamma;

         // Reciprocal lattice vectors
         kBasis_[0][0] = twoPi / a;
         kBasis_[0][1] = -twoPi * cosGamma / (a * sinGamma);
         kBasis_[0][2] = sinTheta * (cosGamma * sinPhi - sinGamma * cosPhi);
         kBasis_[0][2] = twoPi * kBasis_[0][2] / (a * sinGamma * cosTheta);

         kBasis_[1][1] = twoPi/(b*sinGamma);
         kBasis_[1][2] = -twoPi*sinPhi*sinTheta/(b*cosTheta*sinGamma);

         kBasis_[2][2] = twoPi/(c*cosTheta);

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
      } else
      if (lattice == UnitCell<3>::Null) {
         out << "Null";
      } else {
         UTIL_THROW("This should never happen");
      }
      return out;
   }

   /*
   * Assignment operator.
   */
   UnitCell<3>& UnitCell<3>::operator = (UnitCell<3> const & other)
   {
      if (lattice_ != UnitCell<3>::Null) {
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
   * Set lattice system of the unit cell (but not the parameters).
   */
   void UnitCell<3>::set(UnitCell<3>::LatticeSystem lattice)
   {
      UTIL_CHECK(lattice != UnitCell<3>::Null);
      if (lattice_ != UnitCell<3>::Null) {
         UTIL_CHECK(lattice == lattice_);
      }
      isInitialized_ = false;
      lattice_ = lattice;
      setNParameter();
   }

   /*
   * Set state of the unit cell. 
   */
   void UnitCell<3>::set(UnitCell<3>::LatticeSystem lattice,
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
   * Get the length of the unit cell. 
   */
   double UnitCell<3>::volume() const
   {
      double v = 0.0;
      v += rBasis_[0][0]*rBasis_[1][1]*rBasis_[2][2];
      v -= rBasis_[0][0]*rBasis_[1][2]*rBasis_[2][1];
      v += rBasis_[0][1]*rBasis_[1][2]*rBasis_[2][0];
      v -= rBasis_[0][1]*rBasis_[1][0]*rBasis_[2][2];
      v += rBasis_[0][2]*rBasis_[1][0]*rBasis_[2][1];
      v -= rBasis_[0][2]*rBasis_[1][1]*rBasis_[2][0];
      return v;
   }

}
}
