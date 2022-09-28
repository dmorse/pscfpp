#ifndef PSCF_UNIT_CELL_TPP
#define PSCF_UNIT_CELL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"
#include <util/format/Dbl.h> 
#include <util/format/Int.h> 
#include <iomanip>

namespace Pscf
{

   using namespace Util;

   template <int D>
   std::istream& operator >> (std::istream& in,
                              UnitCell<D>& cell)
   {
      in >> cell.lattice_;
      cell.setNParameter();
      for (int i = 0; i < cell.nParameter_; ++i) {
         in >> cell.parameters_[i];
      }
      cell.setLattice();
      return in;
   }

   template <int D>
   std::ostream& operator << (std::ostream& out,
                              UnitCell<D> const & cell)
   {
      out << cell.lattice_;
      for (int i = 0; i < cell.nParameter_; ++i) {
         out << Dbl(cell.parameters_[i], 18, 10);
      }
      return out;
   }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive, int D>
   void serialize(Archive& ar, UnitCell<D>& cell, 
                  const unsigned int version)
   {
      serializeEnum(ar, cell.lattice_, version);
      ar & cell.nParameter_;
      for (int i = 0; i < cell.nParameter_; ++i) {
         ar & cell.parameters_[i];
      }
   }

   template <int D>
   void readUnitCellHeader(std::istream& in, UnitCell<D>& cell)
   {
      std::string label;
      in >> label;
      UTIL_CHECK(label == "crystal_system");
      if (cell.lattice_ == UnitCell<D>::LatticeSystem::Null) {
         in >> cell.lattice_;
         cell.setNParameter();
      } else {
         typename UnitCell<D>::LatticeSystem lattice;
         in >> lattice;
         UTIL_CHECK(lattice == cell.lattice_);
      }

      in >> label;
      UTIL_CHECK(label == "N_cell_param");
      int nParameter;
      in >> nParameter;
      UTIL_CHECK(nParameter == cell.nParameter_);

      in >> label;
      UTIL_CHECK(label == "cell_param");
      for (int i = 0; i < cell.nParameter_; ++i) {
         in >> std::setprecision(15) >> cell.parameters_[i];
      }   
      cell.setLattice();
   }

   template <int D>
   void writeUnitCellHeader(std::ostream& out, UnitCell<D> const& cell)
   {
      out << "crystal_system" << std::endl 
          << "              " << cell.lattice_<< std::endl;
      out << "N_cell_param"   <<  std::endl 
          << "              " << cell.nParameter_<< std::endl;
      out << "cell_param    " <<  std::endl;
      for (int i = 0; i < cell.nParameter_; ++i) {
         out << "    " << Dbl(cell.parameters_[i], 18, 10);
      }
      out << std::endl;
   }

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

      in >> label;
      UTIL_CHECK(label == "format");
      in >> ver1 >> ver2;
 
      in >> label;
      UTIL_CHECK(label == "dim");
      int dim;
      in >> dim;
      UTIL_CHECK(dim == D);

      readUnitCellHeader(in, cell);

      in >> label;
      UTIL_CHECK(label == "group_name");
      in >> groupName;

      in >> label;
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
      out << "format " << Int(ver1,3) << " " << Int(ver2,3) <<  std::endl;
      out << "dim" <<  std::endl 
          << "          " << D << std::endl;
      writeUnitCellHeader(out, cell); 
      out << "group_name" << std::endl 
          << "          " << groupName <<  std::endl;
      out << "N_monomer"  << std::endl 
          << "          " << nMonomer << std::endl;
   }

}
#endif
