#ifndef PRDC_UNIT_CELL_TPP
#define PRDC_UNIT_CELL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"

#include <util/format/Dbl.h> 
#include <util/format/Int.h> 

#include <iomanip>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Read UnitCell<D> from a field file header (fortran PSCF format).
   *
   * If the unit cell has a non-null lattice system on entry, the
   * value read from file must match this existing value, or this
   * function throws an exception. If the lattice system is null on
   * entry, the lattice system value is read from file. In either case,
   * unit cell parameters (dimensions and angles) are updated using
   * values read from file.
   *
   * \param  in  input stream
   * \param  cell  UnitCell<D> to be read
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   void readUnitCellHeader(std::istream& in, UnitCell<D>& cell);

   /**
   * Write UnitCell<D> to a field file header (fortran PSCF format).
   *
   * \param out  output stream
   * \param  cell  UnitCell<D> to be written
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   void writeUnitCellHeader(std::ostream& out, UnitCell<D> const& cell);

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
   void writeUnitCellHeader(std::ostream& out, UnitCell<D> const & cell)
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
