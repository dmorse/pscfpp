/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "replicateUnitCell.h"
#include <prdc/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Create a replicated UnitCell<1>.
   */
   template<>
   void replicateUnitCell(IntVec<1> const & replicas,
                          UnitCell<1> const & cellIn,
                          UnitCell<1> & cellOut)
   {
      UTIL_CHECK(cellIn.lattice() != UnitCell<1>::Null);
      FSArray<double, 6> const & paramIn = cellIn.parameters();

      // Choose lattice type and parameter
      UnitCell<1>::LatticeSystem lattice = UnitCell<1>::Null;
      FArray<double, 6> param;
      if (cellIn.lattice() == UnitCell<1>::Lamellar) {
         lattice = UnitCell<1>::Lamellar;
         param[0] = replicas[0]*paramIn[0];
      } else {
         UTIL_THROW("Invalid lattice system value");
      }

      cellOut.set(lattice);
      FSArray<double, 6> paramOut;
      for (int i=0; i < cellOut.nParameter(); ++i) {
         paramOut.append(param[i]);
      }
      cellOut.set(lattice, paramOut);
   }

   /*
   * Create a replicated UnitCell<2>.
   */
   template<>
   void replicateUnitCell(IntVec<2> const & replicas,
                          UnitCell<2> const & cellIn,
                          UnitCell<2> & cellOut)
   {
      UTIL_CHECK(cellIn.lattice() != UnitCell<2>::Null);
      FSArray<double, 6> const & paramIn = cellIn.parameters();

      // Choose lattice type and parameters
      UnitCell<2>::LatticeSystem lattice = UnitCell<2>::Null;
      FArray<double, 6> param;
      if (cellIn.lattice() == UnitCell<2>::Square) {
         UTIL_CHECK(cellIn.nParameter() == 1);
         if (replicas[0] == replicas[1]) {
            lattice = UnitCell<2>::Square;
            param[0] = replicas[0]*paramIn[0];
         } else {
            lattice = UnitCell<2>::Rectangular;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[1]*paramIn[0];
         }
      } else
      if (cellIn.lattice() == UnitCell<2>::Rectangular) {
         UTIL_CHECK(cellIn.nParameter() == 2);
         lattice = UnitCell<2>::Rectangular;
         param[0] = replicas[0]*paramIn[0];
         param[1] = replicas[1]*paramIn[1];
      } else
      if (cellIn.lattice() == UnitCell<2>::Hexagonal) {
         UTIL_CHECK(cellIn.nParameter() == 1);
         if (replicas[0] == replicas[1]) {
            lattice = UnitCell<2>::Hexagonal;
            param[0] = replicas[0]*paramIn[0];
         } else {
            UTIL_THROW("Prohibited unit cell replication");
         }
      } else
      if (cellIn.lattice() == UnitCell<2>::Rhombic) {
         UTIL_CHECK(cellIn.nParameter() == 2);
         if (replicas[0] == replicas[1]) {
            lattice = UnitCell<2>::Rhombic;
            param[0] = replicas[0]*paramIn[0];
            param[1] = paramIn[1];
         } else {
            lattice = UnitCell<2>::Oblique;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[1]*paramIn[0];
            param[2] = paramIn[2];
         }
      } else
      if (cellIn.lattice() == UnitCell<2>::Oblique) {
         UTIL_CHECK(cellIn.nParameter() == 3);
         lattice = UnitCell<2>::Oblique;
         param[0] = replicas[0]*paramIn[0];
         param[1] = replicas[1]*paramIn[1];
         param[2] = paramIn[2];
      } else {
         UTIL_THROW("Invalid lattice system value");
      }

      // Set cellOut object
      cellOut.set(lattice);
      FSArray<double, 6> paramOut;
      for (int i=0; i < cellOut.nParameter(); ++i) {
         paramOut.append(param[i]);
      }
      cellOut.set(lattice, paramOut);
   }

   /*
   * Create a replicated UnitCell<3>.
   */
   template<>
   void replicateUnitCell(IntVec<3> const & replicas,
                          UnitCell<3> const & cellIn,
                          UnitCell<3> & cellOut)
   {
      UTIL_CHECK(cellIn.lattice() != UnitCell<3>::Null);
      FSArray<double, 6> const & paramIn = cellIn.parameters();

      // Choose lattice type and parameters
      UnitCell<3>::LatticeSystem lattice = UnitCell<3>::Null;
      FArray<double, 6> param;
      int nParameter;
      if (cellIn.lattice() == UnitCell<3>::Cubic) {
         UTIL_CHECK(cellIn.nParameter() == 1);
         if (replicas[0] == replicas[1]) {
            if (replicas[1] == replicas[2]) {
               lattice = UnitCell<3>::Cubic;
               param[0] = replicas[0]*paramIn[0];
               nParameter = 1;
            } else {
               lattice = UnitCell<3>::Tetragonal;
               param[0] = replicas[0]*paramIn[0];
               param[1] = replicas[2]*paramIn[0];
               nParameter = 2;
            }
         } else {
            lattice = UnitCell<3>::Orthorhombic;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[1]*paramIn[0];
            param[2] = replicas[2]*paramIn[0];
            nParameter = 3;
         }
      } else
      if (cellIn.lattice() == UnitCell<3>::Tetragonal) {
         UTIL_CHECK(cellIn.nParameter() == 2);
         if (replicas[0] == replicas[1]) {
            lattice = UnitCell<3>::Tetragonal;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[2]*paramIn[1];
            nParameter = 2;
         } else {
            lattice = UnitCell<3>::Orthorhombic;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[1]*paramIn[0];
            param[2] = replicas[2]*paramIn[1];
            nParameter = 3;
         }
      } else
      if (cellIn.lattice() == UnitCell<3>::Orthorhombic) {
         UTIL_CHECK(cellIn.nParameter() == 3);
         lattice = UnitCell<3>::Orthorhombic;
         param[0] = replicas[0]*paramIn[0];
         param[1] = replicas[1]*paramIn[1];
         param[2] = replicas[2]*paramIn[2];
         nParameter = 3;
      } else
      if (cellIn.lattice() == UnitCell<3>::Monoclinic) {
         UTIL_CHECK(cellIn.nParameter() == 4);
         if (replicas[0] == replicas[1]) {
            lattice = UnitCell<3>::Monoclinic;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[1]*paramIn[1];
            param[2] = replicas[2]*paramIn[2];
            param[3] = paramIn[3];
            nParameter = 4;
         } else {
            UTIL_THROW("Prohibited unit cell replication");
         }
      } else
      if (cellIn.lattice() == UnitCell<3>::Triclinic) {
         UTIL_CHECK(cellIn.nParameter() == 6);
         lattice = UnitCell<3>::Triclinic;
         param[0] = replicas[0]*paramIn[0];
         param[1] = replicas[1]*paramIn[1];
         param[2] = replicas[2]*paramIn[2];
         param[3] = paramIn[3];
         param[4] = paramIn[4];
         param[5] = paramIn[5];
         nParameter = 6;
      } else
      if (cellIn.lattice() == UnitCell<3>::Rhombohedral) {
         UTIL_CHECK(cellIn.nParameter() == 2);
         if (replicas[0] == replicas[1] && replicas[1] == replicas[2]) {
            lattice = UnitCell<3>::Triclinic;
            param[0] = replicas[0]*paramIn[0];
            param[1] = paramIn[1];
            nParameter = 2;
         } else {
            UTIL_THROW("Prohibited unit cell replication");
         }
      } else
      if (cellIn.lattice() == UnitCell<3>::Hexagonal) {
         UTIL_CHECK(cellIn.nParameter() == 2);
         if (replicas[0] == replicas[1]) {
            lattice = UnitCell<3>::Hexagonal;
            param[0] = replicas[0]*paramIn[0];
            param[1] = replicas[2]*paramIn[1];
            nParameter = 2;
         } else {
            UTIL_THROW("Prohibited unit cell replication");
         }
      } else {
         UTIL_THROW("Invalid value of cellIn.lattice()");
      }

      // Set cellOut object
      cellOut.set(lattice);
      UTIL_CHECK(cellOut.nParameter() == nParameter);
      FSArray<double, 6> paramOut;
      for (int i=0; i < nParameter; ++i) {
         paramOut.append(param[i]);
      }
      cellOut.set(lattice, paramOut);
   }

}
}
