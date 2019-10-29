#ifndef PSCF_UNIT_CELL_TPP
#define PSCF_UNIT_CELL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"
#include <util/format/Dbl.h> 
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
         in >> std::setprecision(15) >> cell.parameters_[i];
      }
      cell.setLattice();
      return in;

      #if 0
      std::string label;
      int nParametersIn;
  
      in >> label;
      UTIL_CHECK(label == "crystal_system");
      in >> cell.lattice_;
      cell.setNParameter();

      in >> label;
      UTIL_CHECK(label == "N_cell_param");
      in >> nParametersIn;
      UTIL_CHECK(nParametersIn == cell.nParameter_);

      in >> label;
      UTIL_CHECK(label == "cell_param");


      for (int i = 0; i < cell.nParameter_; ++i) {
         in >> cell.parameters_[i];
      }   
      cell.setLattice();
      return in; 
      #endif

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

      #if 0
      out << "crystal_system" <<  std::endl 
          << "              " << cell.lattice_ << std::endl;
      out << "N_cell_param" <<  std::endl 
          << "                   "<< cell.nParameter_<< std::endl;
      out << "cell_param    " <<  std::endl;
      for (int i = 0; i < cell.nParameter_; ++i) {
         out <<"    "<< Dbl(cell.parameters_[i], 18, 10);
      }
      out<<std::endl;
      return out;
      #endif
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

}
#endif
