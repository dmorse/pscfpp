/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpaceSymmetry.tpp"

namespace Pscf
{

   using namespace Util;

   // Explicit specializations D = 1
   
   template<>
   int SpaceSymmetry<1>::determinant() const
   { return R_(0, 0); }

   template<>
   SpaceSymmetry<1>::Rotation SpaceSymmetry<1>::inverseRotation() const
   {
      // Compute and check the determinant
      int det = R_(0,0);
      if (1 != det*det) {
         std::cout << "Determinant (1D symmetry) =" << det << std::endl;
         UTIL_THROW("Invalid SpaceSymmetry<1>: |det| != 1");
      }

      // Compute the inverse matrix
      SpaceSymmetry<1>::Rotation A;
      A(0,0) = R_(0,0);

      return A;
   }

   // Explicit specializations D = 2
   
   template<>
   int SpaceSymmetry<2>::determinant() const
   {  return (R_(0,0)*R(1,1) - R(0,1)*R(1,0)); }
 
   template<>
   SpaceSymmetry<2>::Rotation SpaceSymmetry<2>::inverseRotation() const
   {
      // Compute and check the determinant
      int det = determinant();
      if (1 != det*det) {
         UTIL_THROW("Invalid SpaceSymmetry<2>: |det| != 1");
      }

      // Compute the inverse matrix
      SpaceSymmetry<2>::Rotation A;
      A(0,0) = R_(1,1)/det;
      A(1,1) = R_(0,0)/det;
      A(0,1) = -R_(0,1)/det;
      A(1,0) = -R_(1,0)/det;

      return A;
   }

   // Explicit specializations D = 3
   
   template<>
   int SpaceSymmetry<3>::determinant() const
   {  
      int det = 0;  
      det += R_(0,0)*(R(1,1)*R(2,2) - R(1,2)*R(2,1));
      det += R_(0,1)*(R(1,2)*R(2,0) - R(1,0)*R(2,2));
      det += R_(0,2)*(R(1,0)*R(2,1) - R(1,1)*R(2,0));
      return det;
   }
 
   template<>
   SpaceSymmetry<3>::Rotation SpaceSymmetry<3>::inverseRotation() const
   {
      // Compute and check the determinant
      int det = determinant();
      if (1 != det*det) {
         UTIL_THROW("Invalid SpaceSymmetry<3>: |det| != 1");
      }

      // Compute the inverse matrix
      SpaceSymmetry<3>::Rotation A;
      A(0,0) = (R_(1,1)*R(2,2) - R(1,2)*R(2,1))/det;
      A(1,0) = (R_(1,2)*R(2,0) - R(1,0)*R(2,2))/det;
      A(2,0) = (R_(1,0)*R(2,1) - R(1,1)*R(2,0))/det;

      A(0,1) = (R_(2,1)*R(0,2) - R(2,2)*R(0,1))/det;
      A(1,1) = (R_(2,2)*R(0,0) - R(2,0)*R(0,2))/det;
      A(2,1) = (R_(2,0)*R(0,1) - R(2,1)*R(0,0))/det;

      A(0,2) = (R_(0,1)*R(1,2) - R(0,2)*R(1,1))/det;
      A(1,2) = (R_(0,2)*R(1,0) - R(0,0)*R(1,2))/det;
      A(2,2) = (R_(0,0)*R(1,1) - R(0,1)*R(1,0))/det;

      return A;
   }

   // Explicit instantiation of required class instances
   
   template class SpaceSymmetry<1>;
   template class SpaceSymmetry<2>;
   template class SpaceSymmetry<3>;

}
