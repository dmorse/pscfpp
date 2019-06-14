#include "SpaceSymmetry.h"

#include <util/format/Int.h>

namespace Pscf
{

   using namespace Util;

   // Define static members
   template<> SpaceSymmetry<2> SpaceSymmetry<2>::identity_;
   template<> bool             SpaceSymmetry<2>::hasIdentity_ = false;

   template<> SpaceSymmetry<3> SpaceSymmetry<3>::identity_;
   template<> bool             SpaceSymmetry<3>::hasIdentity_ = false;

   // Explicit instantiation.
   template class SpaceSymmetry<2>;
   template class SpaceSymmetry<3>;

   template<>
   SpaceSymmetry<1>::Rotation SpaceSymmetry<1>::inverseRotation() const
   {
      // Compute and check the determinant
      int det = R_(0,0);
      if (1 != det*det) {
         UTIL_THROW("Invalid SpaceSymmetry<1>: |det| != 1");
      }

      // Compute the inverse matrix
      SpaceSymmetry<1>::Rotation A;
      A(0,0) = R_(1,1);

      return A; 
   }

   template<>
   SpaceSymmetry<2>::Rotation SpaceSymmetry<2>::inverseRotation() const
   {
      // Compute and check the determinant
      int det;  // Determinant
      det = R_(0,0)*R(1,1) - R(0,1)*R(1,0);
      if (1 != det*det) {
         UTIL_THROW("Invalid SpaceSymmetry<2>: |det| != 1");
      }

      // Compute the inverse matrix
      SpaceSymmetry<2>::Rotation A;
      A(0,0) = R_(1,1)/det;
      A(1,1) = R_(0,0)/det;
      A(0,1) = -R_(1,0)/det;
      A(1,0) = -R_(0,1)/det;

      return A;
   }

}
