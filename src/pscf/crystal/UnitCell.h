#ifndef PSCF_UNIT_CELL_H
#define PSCF_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCellBase.h"
#include <util/format/Dbl.h>
#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Base template for UnitCell<D> classes, D=1, 2 or 3.
   *
   * Explicit specializations are provided for D=1, 2, and 3. In
   * each case, class UnitCell<D> is derived from UnitCellBase<D>,
   * and defines an enumeration named LatticeSystem of the types
   * of Bravais lattice systems in D-dimensional space.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D> 
   class UnitCell  : public UnitCellBase<D>
   {};

   /**
   * istream extractor for a UnitCell<D>.
   *
   * \param  in  input stream
   * \param  cell  UnitCell<D> to be read
   * \return modified input stream
   */
   template <int D>
   std::istream& operator >> (std::istream& in,
                              UnitCell<D>& cell);

   /**
   * ostream inserter for a UnitCell<D>::LatticeSystem.
   *
   * \param out  output stream
   * \param  cell  UnitCell<D> to be written
   * \return modified output stream
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, UnitCell<D>& cell);

   /**
   * Serialize to/from an archive.
   *
   * \param ar       archive
   * \param version  archive version id
   */
   template <class Archive, int D>
   void serialize(Archive& ar, UnitCell<D>& cell, const unsigned int version);

   /**
   * 1D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<1> : public UnitCellBase<1>
   {
   public:

      enum LatticeSystem {Lamellar};

   private:

      LatticeSystem lattice_;

      void setNParameter();
      void setLattice();

      template <int D> 
      friend std::ostream& operator << (std::ostream&, UnitCell<D>& );

      template <int D>
      friend std::istream& operator >> (std::istream&, UnitCell<D>& );

      template <class Archive, int D>
      friend void serialize(Archive& , UnitCell<D>& , const unsigned int );

   };

   /**
   * 2D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<2> : public UnitCellBase<2>
   {
   public:

      enum LatticeSystem {Square, Rectangular, Rhombic, Hexagonal, Oblique};

   private:

      LatticeSystem lattice_;

      void setNParameter();
      void setLattice();

      template <int D> 
      friend std::ostream& operator << (std::ostream&, UnitCell<D>& );

      template <int D>
      friend std::istream& operator >> (std::istream&, UnitCell<D>& );

      template <class Archive, int D>
      friend void serialize(Archive& , UnitCell<D>& , const unsigned int );

   };

   /**
   * 3D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<3> : public UnitCellBase<3>
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

   private:

      LatticeSystem lattice_;

      void setNParameter();
      void setLattice();

      template <int D> 
      friend std::ostream& operator << (std::ostream&, UnitCell<D>& );

      template <int D>
      friend std::istream& operator >> (std::istream&, UnitCell<D>& );

      template <class Archive, int D>
      friend void serialize(Archive& , UnitCell<D>& , const unsigned int );

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

   // Unit Cell inserter (>>) and extractor (<<) operators

   // Implementation Template

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
                              UnitCell<D>& cell)
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
   void serialize(Archive& ar, UnitCell<D>& cell, const unsigned int version)
   {
      serializeEnum(ar, cell.lattice_, version);
      ar & cell.nParameter_;
      for (int i = 0; i < cell.nParameter_; ++i) {
         ar & cell.parameters_[i];
      }
   }

}
#endif
