/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"
#include <util/math/Constants.h>

namespace Pscf {
namespace Rpc
{

   // Explicit instantiations of expandRGridDimension for each D

   template <>
   void
   FieldIo<1>::expandRGridDimension(std::ostream &out,
                                    DArray<RField<1> > const & fields,
                                    UnitCell<1> const & unitCell,
                                    int d,
                                    DArray<int> newGridDimensions) const

   {
      // Check validity of expanded dimension d and newGridDimensions
      UTIL_CHECK(d > 1);
      UTIL_CHECK(d <= 3);
      UTIL_CHECK(newGridDimensions.capacity() == (d - 1));

      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Obtain initial dimension of fields
      IntVec<1> meshDimensions = fields[0].meshDimensions();

      // Set up necessary objects
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      FSArray<double, 6> cellParameters;
      cellParameters.append(unitCell.parameter(0));

      if (d == 2) {
         // 1D expanded to 2D

         // Create necessary objects
         DArray<RField<2> > outFields;
         UnitCell<2> cell;
         IntVec<2> dimensions;

         // Set dimensions
         dimensions[0] = meshDimensions[0];
         dimensions[1] = newGridDimensions[0];

         // Assign unit cell
         if (dimensions[0] == dimensions[1]) {
            cell.set(UnitCell<2>::Square, cellParameters);
         } else {
            cellParameters.append((double)dimensions[1]/dimensions[0] 
                                                  * cellParameters[0]);
            cell.set(UnitCell<2>::Rectangular, cellParameters);
         }

         // Allocate and populate outFields
         outFields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            outFields[i].allocate(dimensions);
            int rank = 0;
            for (int j = 0; j < dimensions[1]; ++j) {
               for (int k = 0; k < dimensions[0]; ++k) {
                  outFields[i][rank] = fields[i][k];
                  rank++;
               }
            }
         }

         // Write Header
         Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
         out << "mesh " <<  std::endl
             << "           " << dimensions << std::endl;

         // Write fields
         MeshIterator<2> itr(dimensions);
         for (itr.begin(); !itr.atEnd(); ++itr) {
            for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(outFields[j][itr.rank()], 18, 15);
            }
            out << std::endl;
         }

      } else if (d == 3) {
         // 1D expanded to 3D

         // Create necessary objects
         DArray<RField<3> > outFields;
         UnitCell<3> cell;
         IntVec<3> dimensions;
         int rank = 0;

         // Set dimensions
         dimensions[0] = meshDimensions[0];
         dimensions[1] = newGridDimensions[0];
         dimensions[2] = newGridDimensions[1];

         // Assign unit cell
         if (dimensions[2] == dimensions[1]) {
            if (dimensions[1] == dimensions[0]) {
               cell.set(UnitCell<3>::Cubic, cellParameters);
            } else {
               cellParameters.append((double)dimensions[1]/dimensions[0] 
                                                     * cellParameters[0]);
               cellParameters.append((double)dimensions[2]/dimensions[0] 
                                                     * cellParameters[0]);
               cell.set(UnitCell<3>::Orthorhombic, cellParameters);
            }
         }
         else {
            if (dimensions[1] == dimensions[0]) {
               cellParameters.append((double)dimensions[2]/dimensions[0] 
                                                     * cellParameters[0]);
               cell.set(UnitCell<3>::Tetragonal, cellParameters);
            } else {
               cellParameters.append((double)dimensions[1]/dimensions[0] 
                                                     * cellParameters[0]);
               cellParameters.append((double)dimensions[2]/dimensions[0] 
                                                     * cellParameters[0]);
               cell.set(UnitCell<3>::Orthorhombic, cellParameters);
            }
         }
         
         // Allocate and populate outFields
         outFields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            outFields[i].allocate(dimensions);
         }
         for (int l = 0; l < dimensions[2]; ++l) {
            for (int k = 0; k < dimensions[1]; ++k) {
               for (int j = 0; j < dimensions[0]; ++j) {
                  for (int i = 0; i < nMonomer; ++i) {
                     outFields[i][rank] = fields[i][j];
                  }
                  ++rank;
               }
            }
         }
         
         // Write Header
         Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
         out << "mesh " <<  std::endl
             << "           " << dimensions << std::endl;
               
         // Write fields
         MeshIterator<3> itr(dimensions);
         for (itr.begin(); !itr.atEnd(); ++itr) {
            for (int j = 0; j < nMonomer; ++j) {
               out << "  " << Dbl(outFields[j][itr.rank()], 18, 15);
            }
            out << std::endl;
         }
      } else {
         UTIL_THROW("Invalid d value");
      }
   }

   template <>
   void
   FieldIo<2>::expandRGridDimension(std::ostream &out,
                                    DArray<RField<2> > const & fields,
                                    UnitCell<2> const & unitCell,
                                    int d,
                                    DArray<int> newGridDimensions) const

   {
      // 2D expanded to 3D

      // Check validity of expanded dimension d and newGridDimensions
      UTIL_CHECK(d == 3);
      UTIL_CHECK(newGridDimensions.capacity() == (d - 2));

      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Obtain initial dimension of fields
      IntVec<2> meshDimensions = fields[0].meshDimensions();

      // Set up necessary objects
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      FSArray<double, 6> cellParameters;
      cellParameters.append(unitCell.parameter(0));

      DArray<RField<3> > outFields;
      UnitCell<3> cell;
      IntVec<3> dimensions;

      // Set dimensions
      dimensions[0] = meshDimensions[0];
      dimensions[1] = meshDimensions[1];
      dimensions[2] = newGridDimensions[0];
      // Set unit cell
      if (unitCell.lattice() == UnitCell<2>::Square) {
         if (newGridDimensions[0] == meshDimensions[0]){
            cell.set(UnitCell<3>::Cubic, cellParameters);
         } else {
            cellParameters.append((double)dimensions[2]/dimensions[0] 
                                                  * cellParameters[0]);
            cell.set(UnitCell<3>::Tetragonal, cellParameters);
         }
      } else if (unitCell.lattice() == UnitCell<2>::Rectangular) {
         cellParameters.append(unitCell.parameter(1));
         cellParameters.append((double)dimensions[2]/dimensions[0] 
                                               * cellParameters[0]);
         cell.set(UnitCell<3>::Orthorhombic, cellParameters);
      } else if (unitCell.lattice() == UnitCell<2>::Hexagonal){
         cellParameters.append((double)dimensions[2]/dimensions[0] 
                                               * cellParameters[0]);
         cell.set(UnitCell<3>::Hexagonal, cellParameters);
      } else if (unitCell.lattice() == UnitCell<2>::Rhombic) {
         cellParameters.append(unitCell.parameter(0));
         cellParameters.append((double)dimensions[2]/dimensions[0] 
                                               * cellParameters[0]);
         cellParameters.append(Constants::Pi / 2);
         cellParameters.append(0.0);
         cellParameters.append(unitCell.parameter(1));
         cell.set(UnitCell<3>::Triclinic, cellParameters);
      } else if (unitCell.lattice() == UnitCell<2>::Oblique) {
         cellParameters.append(unitCell.parameter(1));
         cellParameters.append((double)dimensions[2]/dimensions[0] 
                                               * cellParameters[0]);
         cellParameters.append(Constants::Pi / 2);
         cellParameters.append(0.0);
         cellParameters.append(unitCell.parameter(2));
         cell.set(UnitCell<3>::Triclinic, cellParameters);
      } else {
         UTIL_THROW("Unrecognized 2D lattice system.");
      }

      // Allocate and populate outFields
      outFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         outFields[i].allocate(dimensions);
      }
      
      int q = 0;
      int r = 0;
      int s = 0;
      int n1 =0;
      int n2 =0;
      int n3 =0;

      while (n3 < dimensions[2]) {
         q = 0;
         n2 = 0;
         while (n2 < dimensions[1]) {
            r = q;
            n1 = 0;
            while (n1 < dimensions[0]) {
               for (int i = 0; i < nMonomer; ++i) {
                  outFields[i][s] = fields[i][r];
               }
               r = r + dimensions[1];
               ++s;
               ++n1;
            }
            ++q;
            ++n2;
         }
         ++n3;
      }

      // Write Header
      Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
      out << "mesh " <<  std::endl
          << "           " << dimensions << std::endl;

      // Write fields
      MeshIterator<3> itr(dimensions);
      for (itr.begin(); !itr.atEnd(); ++itr) {
         for (int j = 0; j < nMonomer; ++j) {
            out << "  " << Dbl(outFields[j][itr.rank()], 18, 15);
         }
         out << std::endl;
      }
   }

   template <>
   void
   FieldIo<3>::expandRGridDimension(std::ostream &out,
                                    DArray<RField<3> > const & fields,
                                    UnitCell<3> const & unitCell,
                                    int d,
                                    DArray<int> newGridDimensions) const

   {  UTIL_THROW("expandRGridDimension is invalid when D = 3."); }

   // Class declarations
   template class FieldIo<1>;
   template class FieldIo<2>;
   template class FieldIo<3>;

} // namespace Pscf::Rpc
} // namespace Pscf
