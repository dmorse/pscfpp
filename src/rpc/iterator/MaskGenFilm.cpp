/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilm.tpp"
#include "Iterator.h"

namespace Pscf {
namespace Rpc {

   // Explicit Specializations for setFlexibleParams

   /*
   * Construct an array indicating whether each lattice parameter is 
   * flexible, based on normalVecId and unitCell definitions in param 
   * file as well as the optional user input flexibleParams. Store this
   * array in flexibleParams_ member of the iterator.
   */
   template <>
   void MaskGenFilm<1>::setFlexibleParams()
   {
      // Create flexibleParams array
      FSArray<bool,6> params; 
      params.append(false); // parameter is not flexible

      if (system().iterator().flexibleParams()[0]) {
         Log::file() 
            << "Warning - The lattice parameter is not allowed "
            << "to be flexible for a 1D thin film system."
            << std::endl << std::endl;
      }

      // Pass params into the iterator member of this object
      system().iterator().setFlexibleParams(params);
   }

   template <>
   void MaskGenFilm<2>::setFlexibleParams() 
   {
      if (system().iterator().nFlexibleParams() == 0) { // already rigid
         return;
      }
      
      FSArray<bool,6> params;
      FSArray<bool,6> current = system().iterator().flexibleParams();
      UTIL_CHECK(current.size() == system().unitCell().nParameter());

      // Initialize params to an array of false values
      for (int i = 0; i < current.size(); i++) {
         params.append(false); 
      }

      // In 2D problems, the length of the basis vector with index 
      // normalVecId is fixed and gamma = 90 degrees. If the lattice 
      // system is square, hexagonal, or rhombic, then no parameters
      // can be flexible under these conditions. If rectangular or
      // oblique, then the length of the basis vector parallel to the
      // film may be flexible.
      UnitCell<2>::LatticeSystem lattice = 
                                 system().domain().unitCell().lattice();
      if ((lattice = UnitCell<2>::Rectangular) || 
                                    (lattice = UnitCell<2>::Oblique)) {
         if (normalVecId() == 0) {
            params[1] = current[1];
         } else { // normalVecId() == 1
            params[0] = current[0];
         }
      }

      // Before calling iterator().setFlexibleParams(), check if the number
      // of flexible lattice parameters has changed during this function,
      // and print a warning if so.
      int nFlexParams = 0;
      for (int i = 0; i < params.size(); i++) {
         if (params[i]) nFlexParams++;
      }
      if (nFlexParams < system().iterator().nFlexibleParams()) {
         Log::file() 
            << "***\n"
            << "Notice - Some lattice parameters will be held constant\n"
            << "to comply with the thin film constraint.\n"
            << "***\n" << std::endl;
      }

      // Pass params into the iterator member of this object
      system().iterator().setFlexibleParams(params);
   }

   template <>
   void MaskGenFilm<3>::setFlexibleParams()
   {
      if (system().iterator().nFlexibleParams() == 0) { // already rigid
         return;
      }
      
      FSArray<bool,6> params;
      FSArray<bool,6> current = system().iterator().flexibleParams();
      UTIL_CHECK(current.size() == system().unitCell().nParameter());

      UnitCell<3>::LatticeSystem lattice = 
                                 system().domain().unitCell().lattice();

      // Initialize params to an array of false values
      for (int i = 0; i < current.size(); i++) {
         params.append(false); 
      }

      // There can be up to 3 flexible lattice parameters in 3D: the length 
      // of the two lattice vectors that are not normalVecId, and the angle 
      // between them. The other two angles must be 90 degrees, and the 
      // length of normalVecId is fixed. The crystal system determines which 
      // parameters are flexible.
      if (lattice == UnitCell<3>::Rhombohedral) {

         Log::file() << "Rhombohedral lattice systems are not compatible "
                     << "with a thin film constraint.\n"
                     << "See thin film documentation for more details.\n";
         UTIL_THROW("Cannot use rhombohedral lattice in a thin film system.");

      } else if (lattice == UnitCell<3>::Hexagonal) {

         UTIL_CHECK(normalVecId() == 2); // this is required for hexagonal
         params[0] = current[0]; // a is the only allowed flexible parameter

      } else if (lattice == UnitCell<3>::Tetragonal) {

         // if normalVecId = 2, the only allowed flexibleParam is 0
         // if normalVecId < 2, the only allowed flexibleParam is 1
         if (normalVecId() == 2) {
            params[0] = current[0];
         } else if (normalVecId() < 2) {
            params[1] = current[1];
         }

      } else if (lattice != UnitCell<3>::Cubic) {
         
         for (int i = 0; i < 3; i++) {
            // This if-statement applies for orthorhombic/monoclinic/triclinic
            if (i != normalVecId()) { 
               params[i] = current[i];
            }
         }

         // beta can be flexible in a monoclinic lattice if normalVecId = 1
         if (lattice == UnitCell<3>::Monoclinic && normalVecId() == 1) {
            params[3] = current[3];
         }

         // The angle between the two lattice basis vectors that are not
         // normalVecId can be flexible in a triclinic lattice
         if (lattice == UnitCell<3>::Triclinic) {
            params[normalVecId()+3] = current[normalVecId()+3];
         }

      }

      // Before updating iterator_.flexibleParams_, check if the number
      // of flexible lattice parameters has changed during this function.
      int nFlexParams = 0;
      for (int i = 0; i < params.size(); i++) {
         if (params[i]) nFlexParams++;
      }
      if (nFlexParams < system().iterator().nFlexibleParams()) {
         Log::file() 
            << "***\n"
            << "Notice - Some lattice parameters will be held constant\n"
            << "to comply with the thin film constraint.\n"
            << "***\n" << std::endl;
      }

      // Pass params into the iterator member of this object
      system().iterator().setFlexibleParams(params);
   }

   // Class declarations
   template class MaskGenFilm<1>;
   template class MaskGenFilm<2>;
   template class MaskGenFilm<3>;
}
}