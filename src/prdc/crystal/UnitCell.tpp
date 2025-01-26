#ifndef PRDC_UNIT_CELL_TPP
#define PRDC_UNIT_CELL_TPP

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

namespace Pscf {
namespace Prdc {

   using namespace Util;

   template <int D>
   std::istream& operator >> (std::istream& in,
                              UnitCell<D>& cell)
   {
      typename UnitCell<D>::LatticeSystem lattice;
      in >> lattice;
      cell.set(lattice);
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
      cell.isInitialized_ = false;
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
}
#endif
