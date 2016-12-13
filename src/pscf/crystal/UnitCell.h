#ifndef PSCF_UNIT_CELL_H
#define PSCF_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCellTmpl.h"
#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Base template for UnitCell<D> classes, D=1, 2 or 3.
   *
   * Explicit specializations are provided for D=1, 2, and 3. In
   * each case, class UnitCell<D> is derived from UnitCellTmpl<D>,
   * and defines an enumeration named LatticeSystem of the types
   * of Bravais lattice systems in D-dimensional space.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D> class UnitCell;

   /**
   * 1D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<1> : public UnitCellTmpl<1>
   {
   public:

      enum LatticeSystem {Lamellar};

      void setNParameter();

   private:

      LatticeSystem lattice_;

      friend std::ostream& operator << (std::ostream&, UnitCell<1>& );
      friend std::istream& operator >> (std::istream&, UnitCell<1>& );
   };

   /**
   * 2D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<2> : public UnitCellTmpl<2>
   {
   public:

      enum LatticeSystem {Square, Rectangular, Rhombic, Hexagonal, Oblique};

      void setNParameter();

   private:

      LatticeSystem lattice_;

      friend std::ostream& operator << (std::ostream&, UnitCell<2>& );
      friend std::istream& operator >> (std::istream&, UnitCell<2>& );
   };

   /**
   * 3D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<3> : public UnitCellTmpl<3>
   {
   public:

      /**
      * Enumeration of the 7 possible Bravais lattice systems.
      *
      * Allowed values are: Cubic, Tetragonal, Orthorhombic,
      * Monoclinic, Triclinic, Rhombohedral, and Hexagonal.
      *
      * \ingroup Crystal_Module
      */
      enum LatticeSystem {Cubic, Tetragonal, Orthorhombic, Monoclinic,
                          Triclinic, Rhombohedral, Hexagonal};

      void setNParameter();

   private:

      LatticeSystem lattice_;

      friend std::ostream& operator << (std::ostream&, UnitCell<3>& );
      friend std::istream& operator >> (std::istream&, UnitCell<3>& );
   };

   // Inserter and Extractor Function Declarations

   /**
   * istream extractor for a 1D UnitCell<1>::LatticeSystem.
   *
   * \param in  input stream
   * \param lattice  UnitCell<1>::LatticeSystem to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<1>::LatticeSystem& lattice);

   /**
   * ostream inserter for a 1D UnitCell<1>::LatticeSystem.
   *
   * \param out  output stream
   * \param lattice  UnitCell<1>::LatticeSystem to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<1>::LatticeSystem lattice);

   /**
   * istream extractor for a 2D UnitCell<2>::LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  UnitCell<2>::LatticeSystem to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<2>::LatticeSystem& lattice);

   /**
   * ostream inserter for a 2D UnitCell<2>::LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  UnitCell<2>::LatticeSystem to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<2>::LatticeSystem lattice);

   /**
   * istream extractor for a 3D UnitCell<3>::LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<3>::LatticeSystem& lattice);

   /**
   * ostream inserter for an 3D UnitCell<3>::LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<3>::LatticeSystem lattice);

   /**
   * ostream inserter for a UnitCell<D>::LatticeSystem.
   *
   * \param out  output stream
   * \param lattice  UnitCell<1>::LatticeSystem to be written
   * \return modified output stream
   */
   template <int D>
   std::ostream& operator << (std::ostream& out,
                              UnitCell<D>& cell);

   /**
   * istream extractor for a UnitCell<D>.
   *
   * \param  in       input stream
   * \param  lattice  UnitCell<2>::LatticeSystem to be read
   * \return modified input stream
   */
   template <int D>
   std::istream& operator >> (std::istream& in,
                              UnitCell<D>& cell);

   // Implementation Template

   template <int D>
   std::ostream& operator << (std::ostream& out,
                              UnitCell<D>& cell)
   {
      out << cell.lattice_ << cell.nParameter_;
      for (int i = 0; i < cell.nParameter_; ++i) {
         out << cell.parameters_[i];
      }
      return out;
   }

   template <int D>
   std::istream& operator >> (std::istream& in,
                              UnitCell<D>& cell)
   {
      cell.setNParameter(in);
      for (int i = 0; i < cell.nParameter_; ++i) {
         in >> cell.parameters_[i];
      }
      return in;
   }

}
#endif
