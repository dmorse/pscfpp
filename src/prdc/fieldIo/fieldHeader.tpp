#ifndef PRDC_FIELD_HEADER_TPP
#define PRDC_FIELD_HEADER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/unitCellHeader.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/arithmetic.h>

#include <util/misc/Log.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/math/Constants.h>

#include <string>

namespace Pscf {
namespace Prdc {

   /*
   * Read common part of field header (fortran PSCF format).
   */
   template <int D>
   void readFieldHeader(std::istream& in, 
                        int& ver1, int& ver2, 
                        UnitCell<D>& cell, 
                        std::string& groupName,
                        int& nMonomer)
   {
      std::string label;

      // Read format v1 v2
      in >> label;
      UTIL_CHECK(label == "format");
      in >> ver1;
      in >> ver2;

      // Read dimension of space 
      in >> label;
      UTIL_CHECK(label == "dim");
      int dim;
      in >> dim;
      UTIL_CHECK(dim == D);

      // Read unit cell lattice type, number of parameters, parameters
      // Function declared in UnitCell.h and defined in UnitCell.tpp
      readUnitCellHeader(in, cell);

      // Optionally read groupName
      in >> label;
      groupName = "";
      if (label == "group_name") {
         in >> groupName;
         in >> label;
      }

      // Read nMonomer value
      UTIL_CHECK(label == "N_monomer");
      in >> nMonomer;
      UTIL_CHECK(nMonomer > 0);
   }

   /*
   * Write common part of field header (fortran PSCF format).
   */
   template <int D>
   void writeFieldHeader(std::ostream &out, 
                         int ver1, int ver2,
                         UnitCell<D> const & cell,
                         std::string const & groupName,
                         int nMonomer)
   {
      // Write file format id and dimension of space D
      out << "format " << Int(ver1,3) << " " << Int(ver2,3) <<  std::endl;
      out << "dim" <<  std::endl 
          << "          " << D << std::endl;

      // Write unit cell lattice type, number of parameters, parameters
      // Function declared in UnitCell.h and defined in UnitCell.tpp
      writeUnitCellHeader(out, cell); 

      // Write group_name if not empty
      if (groupName != "") {
         out << "group_name" << std::endl 
             << "          " << groupName <<  std::endl;
      }
  
      // Write nMonomer, number of monomer types
      if (nMonomer != 0) {
         UTIL_CHECK(nMonomer > 0);
         out << "N_monomer"  << std::endl 
             << "          " << nMonomer << std::endl;
      }

      // Note: The option to not write nMonomer when the value is zero
      // is designed to allow this function to be used in contexts in 
      // which the value of N_monomer is irrelevant and possibly unknown,
      // such as files that contain information about stars in a basis 
      // or group operations.

   }

   template <int D>
   void readMeshDimensions(std::istream& in,
                           IntVec<D> const& meshDimensions) 
   {
      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      UTIL_ASSERT(in.good());
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> meshDimensionsIn;
      in >> meshDimensionsIn;
      UTIL_CHECK(in.good());
      if (meshDimensionsIn != meshDimensions) {
         Log::file()
           << "Inconsistent mesh dimensions:\n"
           << "meshDimensions (expected)  = " << meshDimensions << "\n"
           << "meshDimensions (from file) = " << meshDimensionsIn << "\n";
         UTIL_THROW("Unexpected mesh dimensions in field file header");
      }
   }

   template <int D>
   void writeMeshDimensions(std::ostream &out,
                            IntVec<D> const& meshDimensions)
   {
      out << "mesh " <<  std::endl
          << "           " << meshDimensions << std::endl;
   }

} // namespace Prdc
} // namespace Pscf
#endif
