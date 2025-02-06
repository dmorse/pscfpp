/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilmBase.tpp"
#include <prdc/crystal/paramIdConversions.h>

namespace Pscf {
namespace Prdc
{

   using namespace Util;

   // Explicit Specializations for checkLatticeVectors
   
   /*
   * Check that user-defined lattice basis vectors are compatible with 
   * the thin film constraint.
   * 
   * In 1D, there is nothing to do; the lattice basis vector is correct in 
   * all cases.
   */
   template <>
   void MaskGenFilmBase<1>::checkLatticeVectors() const 
   {} // do nothing

   
   // In 2D, we require that gamma = 90°.
   template <>
   void MaskGenFilmBase<2>::checkLatticeVectors() const 
   {
      RealVec<2> a, b;
      a = systemLatticeVector(0);
      b = systemLatticeVector(1);

      double gamma = dot(a,b);
      if (gamma > 1e-8) { // Dot product between a and b should be 0
         UTIL_THROW("ERROR: Lattice basis vectors must be orthogonal");
      }
   } 

   /*
   * In 3D, we require that there be one lattice basis vector that is 
   * orthogonal to the walls (parameter with index normalVecId), and two
   * that are parallel to the walls.
   */
   template <>
   void MaskGenFilmBase<3>::checkLatticeVectors() const 
   {
      RealVec<3> a, b, c;
      a = systemLatticeVector(0);
      b = systemLatticeVector(1);
      c = systemLatticeVector(2);

      double alpha, beta, gamma;
      gamma = dot(a,b);
      beta = dot(a,c);
      alpha = dot(b,c);

      if (normalVecId() == 0) {
         if (beta > 1e-8 || gamma > 1e-8) {
            UTIL_THROW("ERROR: beta and gamma must be 90 degrees");
         }
      } else if (normalVecId() == 1) {
         if (alpha > 1e-8 || gamma > 1e-8) {
            UTIL_THROW("ERROR: alpha and gamma must be 90 degrees");
         }
      } else { // normalVecId == 2
         if (alpha > 1e-8 || beta > 1e-8) {
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
   MaskGenFilmBase<1>::modifyFlexibleParams(FSArray<bool,6> current,
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
   MaskGenFilmBase<2>::modifyFlexibleParams(FSArray<bool,6> current,
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
   MaskGenFilmBase<3>::modifyFlexibleParams(FSArray<bool,6> current,
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
         // not normalVecId to be flexible
         updated[normalVecId_+3] = current[normalVecId_+3];

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
   template class MaskGenFilmBase<1>;
   template class MaskGenFilmBase<2>;
   template class MaskGenFilmBase<3>;
}
}