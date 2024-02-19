#ifndef RPC_FILM_ITERATOR_TPP
#define RPC_FILM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIterator.h"
#include "rpc/System.h"
#include "prdc/crystal/UnitCell.h"
#include "pscf/math/RealVec.h"

namespace Pscf {
namespace Rpc
{
   using namespace Util;
   using namespace Pscf::Prdc;

   // 1D Partial Specializations

   /*
   * Constructor.
   */
   template <typename IteratorType>
   FilmIterator<1, IteratorType>::FilmIterator(System<1>& system)
    : FilmIteratorBase<1, IteratorType>(system)
   {  setClassName(iterator().className().append("Film").c_str()); }

   /*
   * Construct an array indicating whether each lattice parameter is 
   * flexible, based on normalVecId and unitCell definitions in param 
   * file as well as the optional user input flexibleParams. Store this
   * array in flexibleParams_ member of this object, as well the 
   * flexibleParams_ member of the iterator within this object.
   * Uses the flexibleParams_ member of the iterator within this object
   * as a starting point.
   */
   template <typename IteratorType>
   void FilmIterator<1, IteratorType>::setFlexibleParams()
   {
      // Create flexibleParams array
      FSArray<bool,6> params; 
      params.append(false); // parameter is not flexible

      if (iterator().flexibleParams()[0]) {
         Log::file() 
            << "Warning - The lattice parameter is not allowed "
            << "to be flexible for a 1D thin film system."
            << std::endl;
      }

      // Store params in flexibleParams_ member of this object
      setFlexibleParams(params);

      // Pass params into the iterator member of this object
      iterator().setFlexibleParams(params);
   }

   /*
   * Check that user-defined lattice basis vectors are compatible with 
   * the thin film constraint (in 1D, there is nothing to do; the lattice
   * basis vector is correct in all cases).
   */
   template <typename IteratorType>
   void FilmIterator<1, IteratorType>::checkLatticeVectors() const 
   {} // do nothing


   // 2D Partial Specializations

   /*
   * Constructor.
   */
   template <typename IteratorType>
   FilmIterator<2, IteratorType>::FilmIterator(System<2>& system)
    : FilmIteratorBase<2, IteratorType>(system)
   {  setClassName(iterator().className().append("Film").c_str()); }

   /*
   * Construct an array indicating whether each lattice parameter is
   * flexible, based on normalVecId and unitCell definitions in param 
   * file as well as the optional user input flexibleParams. Store this
   * array in flexibleParams_ member of this object, as well the 
   * flexibleParams_ member of the iterator within this object.
   * Uses the flexibleParams_ member of the iterator within this object
   * as a starting point.
   */
   template <typename IteratorType>
   void FilmIterator<2, IteratorType>::setFlexibleParams() 
   {
      FSArray<bool,6> params;
      FSArray<bool,6> current = iterator().flexibleParams();
      UTIL_CHECK(current.size() == system().unitCell().nParameter());

      if (iterator().nFlexibleParams() == 0) { // lattice is already rigid
         setFlexibleParams(current);
         return;
      }

      // Initialize params to an array of false values
      for (int i = 0; i < current.size(); i++) {
         params.append(false); 
      }

      // If the lattice system is square, no params are flexible. Otherwise,
      // params should, at most, contain only the index of the lattice basis
      // vector that is not normalVecId. The length of normalVecId is fixed
      //and gamma = 90 degrees for all 2D thin film calculations.
      if (system().domain().unitCell().lattice() != UnitCell<2>::Square) {
         for (int i = 0; i < params.size(); i++) {
            if ((i != normalVecId()) && (i < 2)) {
               params[i] = true;
            }
         }
      }

      // Store params in flexibleParams_ member of this object
      setFlexibleParams(params);

      // Before updating iterator_.flexibleParams_, check if the number
      // of flexible lattice parameters has changed during this function.
      if (nFlexibleParams() < iterator().nFlexibleParams()) {
         Log::file() 
            << "***Notice - Some lattice parameters will be held constant\n"
            << "to comply with the thin film constraint.***"
            << std::endl;
      }

      // Pass params into the iterator member of this object
      iterator().setFlexibleParams(params);
   }

   /*
   * Check that user-defined lattice basis vectors are compatible with 
   * the thin film constraint (in 2D, we require that gamma = 90°)
   */
   template <typename IteratorType>
   void FilmIterator<2, IteratorType>::checkLatticeVectors() const 
   {
      RealVec<2> a, b;
      a = system().domain().unitCell().rBasis(0);
      b = system().domain().unitCell().rBasis(1);

      double gamma = dot(a,b);
      if (gamma > 1e-8) { // Dot product between a and b should be 0
         UTIL_THROW("ERROR: Lattice basis vectors must be orthogonal when wall is present");
      }
   } 


   // 3D Partial Specializations

   /*
   * Constructor.
   */
   template <typename IteratorType>
   FilmIterator<3, IteratorType>::FilmIterator(System<3>& system)
    : FilmIteratorBase<3, IteratorType>(system)
   {  setClassName(iterator().className().append("Film").c_str()); }

   /*
   * Construct an array indicating whether each lattice parameter is 
   * flexible, based on normalVecId and unitCell definitions in param 
   * file as well as the optional user input flexibleParams. Store this
   * array in flexibleParams_ member of this object, as well the 
   * flexibleParams_ member of the iterator within this object.
   * Uses the flexibleParams_ member of the iterator within this object
   * as a starting point.
   */
   template <typename IteratorType>
   void FilmIterator<3, IteratorType>::setFlexibleParams()
   {
      FSArray<bool,6> params;
      FSArray<bool,6> current = iterator().flexibleParams();
      UTIL_CHECK(current.size() == system().unitCell().nParameter());

      UnitCell<3>::LatticeSystem lattice = 
                                 system().domain().unitCell().lattice();

      // Initialize params to an array of false values
      for (int i = 0; i < current.size(); i++) {
         params.append(false); 
      }

      if (iterator().nFlexibleParams() == 0) { // lattice is already rigid
         setFlexibleParams(params);
         return;
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
         if (current[0]) {
            params[0] = true; // a is the only allowed flexible parameter
         }

      } else if (lattice == UnitCell<3>::Tetragonal) {

         // if normalVecId = 2, the only allowed flexibleParam is 0
         // if normalVecId < 2, the only allowed flexibleParam is 1
         if (current[0] && normalVecId() == 2) {
            params[0] = true;
         } else if (current[1] && normalVecId() < 2) {
            params[1] = true;
         }

      } else if (lattice != UnitCell<3>::Cubic) {
         
         for (int i = 0; i < 3; i++) {
            // This if-statement applies for orthorhombic/monoclinic/triclinic
            if (i != normalVecId() && current[i]) { 
               params[i] = true;
            }
         }

         // beta can be flexible in a monoclinic lattice if normalVecId = 1
         if (lattice == UnitCell<3>::Monoclinic && normalVecId() == 1
               && current[3]) {
            params[3] = true;
         }

         // The angle between the two lattice basis vectors that are not
         // normalVecId can be flexible in a triclinic lattice
         if (lattice == UnitCell<3>::Triclinic && current[normalVecId()+3]) {
            params[normalVecId()+3] = true;
         }

      }

      // Store params in flexibleParams_ member of this object
      setFlexibleParams(params);

      // Before updating iterator_.flexibleParams_, check if the number
      // of flexible lattice parameters has changed during this function.
      if (nFlexibleParams() < iterator().nFlexibleParams()) {
         Log::file() 
            << "***Notice - Some lattice parameters will be held constant\n"
            << "to comply with the thin film constraint.***"
            << std::endl;
      }

      // Pass params into the iterator member of this object
      iterator().setFlexibleParams(params);
   }

   /*
   * Check that user-defined lattice basis vectors are compatible with 
   * the thin film constraint (in 3D, we require that there be one 
   * lattice basis vector that is orthogonal to the walls, and two that
   * are parallel to the walls; the orthogonal vector is normalVecId).
   */
   template <typename IteratorType>
   void FilmIterator<3, IteratorType>::checkLatticeVectors() const 
   {
      RealVec<3> a, b, c;
      a = system().domain().unitCell().rBasis(0);
      b = system().domain().unitCell().rBasis(1);
      c = system().domain().unitCell().rBasis(2);
      double alpha, beta, gamma;
      gamma = dot(a,b);
      beta = dot(a,c);
      alpha = dot(b,c);

      if (normalVecId() == 0) {
         if (beta > 1e-8 || gamma > 1e-8) {
            UTIL_THROW("ERROR: If normalVecId = 0, beta and gamma must be 90 degrees");
         }
      } else if (normalVecId() == 1) {
         if (alpha > 1e-8 || gamma > 1e-8) {
            UTIL_THROW("ERROR: If normalVecId = 1, alpha and gamma must be 90 degrees");
         }
      } else { // normalVecId == 2
         if (alpha > 1e-8 || beta > 1e-8) {
            UTIL_THROW("ERROR: If normalVecId = 2, alpha and beta must be 90 degrees");
         }
      }
   }

} // namespace Rpc
} // namespace Pscf
#endif
