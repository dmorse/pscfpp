/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenMaskBase.tpp"
#include <pscf/crystal/paramIdConversions.h>

namespace Pscf {
namespace Rp
{

   using namespace Util;

   namespace {

      /*
      * Are two lattice basis vectors orthogonal? The test is on the
      * cosine of their angle, so it does not depend on the unit of
      * length (a right angle with 11 digits has |cos| <= 5.0E-11).
      */
      template <int D>
      bool isOrthogonal(RealVec<D> const & a, RealVec<D> const & b)
      {
         double ab = dot(a, b);
         return (std::abs(ab) <= 1.0E-8*std::sqrt(dot(a, a)*dot(b, b)));
      }

   }

   // Explicit Specializations for checkLatticeVectors
   
   /*
   * Check that user-defined lattice basis vectors are compatible with 
   * the thin film constraint.
   * 
   * In 1D, there is nothing to do; the lattice basis vector is correct in 
   * all cases.
   */
   template <>
   void FilmFieldGenMaskBase<1>::checkLatticeVectors() const 
   {} // do nothing

   
   // In 2D, we require that gamma = 90°.
   template <>
   void FilmFieldGenMaskBase<2>::checkLatticeVectors() const 
   {
      RealVec<2> a, b;
      a = systemLatticeVector(0);
      b = systemLatticeVector(1);

      if (!isOrthogonal(a, b)) {
         UTIL_THROW("ERROR: Lattice basis vectors must be orthogonal");
      }
   } 

   /*
   * In 3D, we require that there be one lattice basis vector that is 
   * orthogonal to the walls (parameter with index normalVecId), and two
   * that are parallel to the walls.
   */
   template <>
   void FilmFieldGenMaskBase<3>::checkLatticeVectors() const 
   {
      RealVec<3> a, b, c;
      a = systemLatticeVector(0);
      b = systemLatticeVector(1);
      c = systemLatticeVector(2);

      // Right angles alpha (b,c), beta (a,c) and gamma (a,b)?
      bool alphaOk = isOrthogonal(b, c);
      bool betaOk = isOrthogonal(a, c);
      bool gammaOk = isOrthogonal(a, b);

      if (normalVecId() == 0) {
         if (!betaOk || !gammaOk) {
            UTIL_THROW("ERROR: beta and gamma must be 90 degrees");
         }
      } else if (normalVecId() == 1) {
         if (!alphaOk || !gammaOk) {
            UTIL_THROW("ERROR: alpha and gamma must be 90 degrees");
         }
      } else { // normalVecId == 2
         if (!alphaOk || !betaOk) {
            UTIL_THROW("ERROR: alpha and beta must be 90 degrees");
         }
      }
   }

   // Explicit Specializations for modifyFlexibleParams

   /*
   * Modifies a flexibleParams array to be compatible with this mask. (1D)
   */
   template <>
   FSArray<bool,6> 
   FilmFieldGenMaskBase<1>::modifyFlexibleParams(FSArray<bool,6> current,
                                            UnitCell<1> const & cell) const
   {
      UTIL_CHECK(current.size() == cell.nParameter());
      UTIL_CHECK(cell.nParameter() == 1);

      if ((current[0]) && (!hasFBulk())) {
         // If user specified that the parameter is flexible but 
         // did not provide a value for fBulk, the parameter will 
         // be changed to rigid
         current[0] = false;
         Log::file() 
            << "Warning - The lattice parameter is not allowed "
            << "to be flexible\n for a 1D thin film system "
            << "unless fBulk is provided."
            << std::endl;
      }

      return current;
   }

   /*
   * Modifies a flexibleParams array to be compatible with this mask. (2D)
   */
   template <>
   FSArray<bool,6> 
   FilmFieldGenMaskBase<2>::modifyFlexibleParams(FSArray<bool,6> current,
                                            UnitCell<2> const & cell) const
   {
      UTIL_CHECK(current.size() == cell.nParameter());
      
      FSArray<bool,6> updated = current;

      // In 2D problems, gamma is fixed at 90 degrees, and the length
      // of the basis vector with index normalVecId is fixed as well,
      // unless the input parameter fBulk is provided, in which case 
      // it may be flexible. If fBulk is not provided and the lattice 
      // system is square, hexagonal, or rhombic, then no parameters
      // can be flexible under these conditions. If rectangular or
      // oblique, then the length of the basis vector parallel to the
      // film may be flexible.

      // First, set gamma to rigid
      if (cell.lattice() == UnitCell<2>::Rhombic) updated[1] = false;
      if (cell.lattice() == UnitCell<2>::Oblique) updated[2] = false;

      // Next, set the parameter corresponding to normalVecId to rigid
      // (unless fBulk was provided)
      if (!hasFBulk()) {
         int reducedId = convertFullParamIdToReduced<2>(normalVecId_, 
                                                        cell.lattice());
         updated[reducedId] = false;
      }

      // Check if the number of flexible lattice parameters has changed 
      // during this function, and print a warning if so.
      bool warn = false;
      for (int i = 0; i < updated.size(); i++) {
         if (updated[i] != current[i]) {
            warn = true;
            break;
         }
      }
      if (warn) {
         Log::file() 
            << "***\n"
            << "Notice - Some lattice parameters will be held constant\n"
            << "to comply with the thin film constraint.\n"
            << "***" << std::endl;
      }

      return updated;
   }

   /*
   * Modifies a flexibleParams array to be compatible with this mask. (3D)
   */
   template <>
   FSArray<bool,6> 
   FilmFieldGenMaskBase<3>::modifyFlexibleParams(FSArray<bool,6> current,
                                            UnitCell<3> const & cell) const
   {
      UTIL_CHECK(current.size() == cell.nParameter());
      
      FSArray<bool,6> updated = current;

      // If fBulk is not provided, then there can be up to 3 flexible 
      // lattice parameters in 3D: the length of the two lattice vectors 
      // that are not normalVecId, and the angle between them. The other 
      // two angles must be 90 degrees, and the length of normalVecId is 
      // fixed. The crystal system determines which parameters are flexible.
      // 
      // If fBulk is provided, then the length of the basis vector with
      // index normalVecId may also be flexible. 

      // First, set angles to rigid if necessary
      if (cell.lattice() == UnitCell<3>::Rhombohedral) {

         Log::file() << "Rhombohedral lattice systems are not compatible "
                     << "with a thin film constraint.\n"
                     << "See thin film documentation for more details.\n";
         UTIL_THROW("Cannot use rhombohedral lattice in a thin film system.");

      } else if (cell.lattice() == UnitCell<3>::Monoclinic) {

         // Beta can be flexible in a monoclinic lattice iff normalVecId = 1
         if (normalVecId_ != 1) {
            updated[3] = false;
         }

      } else if (cell.lattice() == UnitCell<3>::Triclinic) {

         // Set all angles to rigid
         updated[3] = false;
         updated[4] = false;
         updated[5] = false;

         // Allow the angle between the two lattice basis vectors that are
         // not normalVecId to be flexible. The triclinic parameters are
         // (a, b, c, phi, theta, gamma): b lies in the x-y plane at angle
         // gamma from a, and c has polar angle theta and azimuth phi.
         //  - normalVecId = 2: theta = 0; the free angle is gamma (5).
         //  - normalVecId = 1: gamma = 90 deg, phi = 0; the angle between
         //    a and c is 90 deg - theta, so the free angle is theta (4).
         //  - normalVecId = 0: gamma = phi = 90 deg; the angle between
         //    b and c is 90 deg - theta, so the free angle is theta (4).
         //    (phi would tilt c out of the plane normal to a.)
         int angleId = (normalVecId_ == 2) ? 5 : 4;
         updated[angleId] = current[angleId];

         // Theta is the in-plane angle only if phi has the value assumed
         // above. At theta = 0 the orthogonality check passes for any
         // phi, but relaxing theta would then tilt c toward the wall
         // normal. Require phi explicitly when theta is flexible.
         if (angleId == 4 && updated[4]) {
            double phi = cell.parameter(3);
            double t = (normalVecId_ == 0) ? std::cos(phi) : std::sin(phi);
            if (std::abs(t) > 1.0E-8) {
               UTIL_THROW("Flexible theta in a triclinic film needs phi = "
                          "pi/2 (normalVecId 0) or 0 (normalVecId 1)");
            }
         }

      }
      
      // Next, set the parameter corresponding to normalVecId to rigid
      // (unless fBulk was provided)
      if (!hasFBulk()) {
         int reducedId = convertFullParamIdToReduced<3>(normalVecId_, 
                                                        cell.lattice());
         updated[reducedId] = false;
      }

      // Check if the number of flexible lattice parameters has changed 
      // during this function, and print a warning if so.
      bool warn = false;
      for (int i = 0; i < updated.size(); i++) {
         if (updated[i] != current[i]) {
            warn = true;
            break;
         }
      }
      if (warn) {
         Log::file() 
            << "***\n"
            << "Notice - Some lattice parameters will be held constant\n"
            << "to comply with the thin film constraint.\n"
            << "***" << std::endl;
      }

      return updated;
   }

   // Class declarations
   template class FilmFieldGenMaskBase<1>;
   template class FilmFieldGenMaskBase<2>;
   template class FilmFieldGenMaskBase<3>;
}
}
