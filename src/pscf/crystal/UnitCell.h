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
#include <iomanip>

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
   class UnitCell : public UnitCellBase<D>
   {};

   // Function declations (friends of explicit instantiations)

   /**
   * istream input extractor for a UnitCell<D>.
   *
   * \param  in  input stream
   * \param  cell  UnitCell<D> to be read
   * \return modified input stream
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   std::istream& operator >> (std::istream& in, UnitCell<D>& cell);

   /**
   * ostream output inserter for a UnitCell<D>.
   *
   * \param out  output stream
   * \param  cell  UnitCell<D> to be written
   * \return modified output stream
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, UnitCell<D> const& cell);

   /**
   * Serialize to/from an archive.
   *
   * \param ar       archive
   * \param version  archive version id
   * \ingroup Pscf_Crystal_Module
   */
   template <class Archive, int D>
   void serialize(Archive& ar, UnitCell<D>& cell,
                  const unsigned int version);

   /**
   * Read UnitCell<D> from a field file header (fortran pscf format).
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
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   void readUnitCellHeader(std::istream& in, UnitCell<D>& cell);

   /**
   * Write UnitCell<D> to a field file header (fortran pscf format).
   *
   * \param out  output stream
   * \param  cell  UnitCell<D> to be written
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   void writeUnitCellHeader(std::ostream& out, UnitCell<D> const& cell);

   // 1D Unit Cell

   /**
   * 1D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<1> : public UnitCellBase<1>
   {
   public:

      /**
      * Enumeration of 1D lattice system types.
      */
      enum LatticeSystem {Lamellar, Null};

      /**
      * Constructor
      */
      UnitCell();

   private:

      LatticeSystem lattice_;

      void setNParameter();

      void setBasis();

   // friends:

      template <int D>
      friend std::istream& operator >> (std::istream&, UnitCell<D>& );

      template <int D>
      friend std::ostream& operator << (std::ostream&, UnitCell<D> const&);

      template <class Archive, int D>
      friend void serialize(Archive& , UnitCell<D>& , const unsigned int);

      template <int D>
      friend void readUnitCellHeader(std::istream&, UnitCell<D>& );

      template <int D>
      friend void writeUnitCellHeader(std::ostream&, UnitCell<D> const&);

   };

   /**
   * istream extractor for a 1D UnitCell<1>::LatticeSystem.
   *
   * \param in  input stream
   * \param lattice  UnitCell<1>::LatticeSystem to be read
   * \return modified input stream
   * \ingroup Pscf_Crystal_Module
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<1>::LatticeSystem& lattice);

   /**
   * ostream inserter for a 1D UnitCell<1>::LatticeSystem.
   *
   * \param out  output stream
   * \param lattice  UnitCell<1>::LatticeSystem to be written
   * \return modified output stream
   * \ingroup Pscf_Crystal_Module
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<1>::LatticeSystem lattice);

   /**
   * 2D crystal unit cell.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <>
   class UnitCell<2> : public UnitCellBase<2>
   {
   public:

      /**
      * Enumeration of 2D lattice system types.
      */
      enum LatticeSystem {Square, Rectangular, Rhombic, Hexagonal,
                          Oblique, Null};

      /**
      * Constructor
      */
      UnitCell();

   private:

      /**
      * Lattice system (square, rectangular, etc.)
      */
      LatticeSystem lattice_;

      void setNParameter();

      void setBasis();

   // friends:

      template <int D>
      friend std::istream& operator >> (std::istream&, UnitCell<D>& );

      template <int D>
      friend std::ostream& operator << (std::ostream&, UnitCell<D> const&);

      template <class Archive, int D>
      friend void serialize(Archive& , UnitCell<D>& , const unsigned int );

      template <int D>
      friend void readUnitCellHeader(std::istream&, UnitCell<D>& );

      template <int D>
      friend void writeUnitCellHeader(std::ostream&, UnitCell<D> const&);

   };

   /**
   * istream extractor for a 2D UnitCell<2>::LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  UnitCell<2>::LatticeSystem to be read
   * \return modified input stream
   * \ingroup Pscf_Crystal_Module
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

   // 3D crystal unit cell

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
      * Enumeration of the 7 possible 3D Bravais lattice systems.
      *
      * Allowed non-null values are: Cubic, Tetragonal, Orthorhombic,
      * Monoclinic, Triclinic, Rhombohedral, and Hexagonal.
      */
      enum LatticeSystem {Cubic, Tetragonal, Orthorhombic, Monoclinic,
                          Triclinic, Rhombohedral, Hexagonal, Null};

      /**
      * Constructor
      */
      UnitCell();

   private:

      LatticeSystem lattice_;

      void setNParameter();

      void setBasis();

   // friends:

      template <int D>
      friend std::istream& operator >> (std::istream&, UnitCell<D>& );

      template <int D>
      friend std::ostream& operator << (std::ostream&, UnitCell<D> const&);

      template <class Archive, int D>
      friend void serialize(Archive& , UnitCell<D>& , const unsigned int);

      template <int D>
      friend void readUnitCellHeader(std::istream&, UnitCell<D>& );

      template <int D>
      friend void writeUnitCellHeader(std::ostream&, UnitCell<D> const&);

   };

   /**
   * istream extractor for a 3D UnitCell<3>::LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be read
   * \return modified input stream
   * \ingroup Pscf_Crystal_Module
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<3>::LatticeSystem& lattice);

   /**
   * ostream inserter for an 3D UnitCell<3>::LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be written
   * \return modified output stream
   * \ingroup Pscf_Crystal_Module
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<3>::LatticeSystem lattice);

}

#include "UnitCell.tpp"
#endif
