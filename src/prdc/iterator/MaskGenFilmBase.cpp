/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilmBase.tpp"

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
      a = getLatticeVector(0);
      b = getLatticeVector(1);

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
      a = getLatticeVector(0);
      b = getLatticeVector(1);
      c = getLatticeVector(2);

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

   // Class declarations
   template class MaskGenFilmBase<1>;
   template class MaskGenFilmBase<2>;
   template class MaskGenFilmBase<3>;
}
}