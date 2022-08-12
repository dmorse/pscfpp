#ifndef PSPC_FILM_ITERATOR_TPP
#define PSPC_FILM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIterator.h"
#include "pscf/crystal/UnitCell.h"
#include "pscf/math/RealVec.h"
#include "pspc/System.h"

namespace Pscf {
namespace Pspc
{
   using namespace Util;

   // 1D Partial Specializations

   /*
   * Constructor.
   */
   template <typename IteratorType>
   FilmIterator<1, IteratorType>::FilmIterator(System<1>& system)
    : FilmIteratorBase<1, IteratorType>(system)
   {  setClassName(iterator().className().append("Film").c_str()); }

   /*
   * Output the indices of each flexible lattice parameter, based on
   * normalVecId and unitCell definitions in param file
   */
   template <typename IteratorType>
   void FilmIterator<1, IteratorType>::setFlexibleParams()
   {
      FSArray<int, 6> params; // empty FSArray, no flexible parameters

      Log::file() << std::endl 
         << "Warning - The lattice parameter is not allowed "
         << "to be flexible for a 1D system that contains walls.\n"
         << std::endl;

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
   * Output the indices of each flexible lattice parameter, based on
   * normalVecId and unitCell definitions in param file
   */
   template <typename IteratorType>
   void FilmIterator<2, IteratorType>::setFlexibleParams() 
   {
      FSArray<int, 6> params;

      // If the lattice system is square, params should be empty. Otherwise,
      // params should contain only the length of the vector that is not 
      // normalVecId. The length of normalVecId is fixed, and gamma = 90 degrees.
      if (system().domain().unitCell().lattice() != UnitCell<2>::Square) {
         if (normalVecId() == 0) {
            params.append(1);
         } else { // normalVecId() == 1
            params.append(0);
         }
      }

      if (params.size() == 0) {
         Log::file() << std::endl 
            << "Warning - None of the lattice parameters are allowed\n"
            << "to be flexible for this choice of lattice system\n"
            << "when the system is confined in a thin film."
            << std::endl;
      }

      // Store params in flexibleParams_ member of this object
      setFlexibleParams(params);

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
   * Output the indices of each flexible lattice parameter, based on
   * normalVecId and unitCell definitions in param file
   */
   template <typename IteratorType>
   void FilmIterator<3, IteratorType>::setFlexibleParams()
   {
      FSArray<int, 6> params;
      UnitCell<3>::LatticeSystem lattice = system().domain().unitCell().lattice();

      // params can contain up to 3 lattice parameters: the length of the two 
      // lattice vectors that are not normalVecId, and the angle between them.
      // The other two angles must be 90 degrees, and the length of normalVecId
      // is fixed. The crystal system determines which parameters are flexible.
      if (lattice == UnitCell<3>::Hexagonal || 
         lattice == UnitCell<3>::Rhombohedral) {
         UTIL_CHECK(normalVecId() == 2); // this is required for hex/rhombohedral
         params.append(0);
      } else if (lattice == UnitCell<3>::Tetragonal) {
         if (normalVecId() == 2) {
            params.append(0);
         } else {
            params.append(1);
         }
      } else if (lattice != UnitCell<3>::Cubic) {
         // The for-loop below applies to orthorhombic, monoclinic, & triclinic
         for (int i = 0; i < 3; i++) {
            if (normalVecId() != i) {
               params.append(i);
            }
         }
         if (lattice == UnitCell<3>::Monoclinic && normalVecId() == 1) {
            params.append(3); // beta is flexible if normalVecId == 1
         } else if (lattice == UnitCell<3>::Triclinic) {
            params.append(normalVecId() + 3);
         }
      }

      if (params.size() == 0) {
         Log::file() << std::endl 
            << "Warning - None of the lattice parameters are allowed\n"
            << "to be flexible for this choice of lattice system\n"
            << "when the system is confined in a thin film."
            << std::endl;
      }

      // Store params in flexibleParams_ member of this object
      setFlexibleParams(params);

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

} // namespace Pspc
} // namespace Pscf
#endif