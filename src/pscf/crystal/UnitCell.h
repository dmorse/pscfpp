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
   * Explicit specializations are provided for D=1, 2, and 3.
   * In each case UnitCell<D> is derived from UnitCellTmpl<D>,
   * and defines an enumeration named LatticeSystem of types
   * of Bravais lattice systems for the appropriate dimension.
   */
   template <int D> class UnitCell;

   /**
   * 1D crystal unit cell.
   *
   * \ingroup Pscf_Base_Module
   */
   template <>
   class UnitCell<1> : public UnitCellTmpl<1>
   {
   public:

      enum LatticeSystem {Lamellar};

   };
   
   /**
   * 2D crystal unit cell.
   *
   * \ingroup Pscf_Base_Module
   */
   template <>
   class UnitCell<2> : public UnitCellTmpl<2>
   {
   public:

      enum LatticeSystem {Square, Rectangular, Rhombic, Hexagonal, Oblique};

   };
   
   /**
   * 3D crystal unit cell.
   *
   * \ingroup Pscf_Base_Module
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

   };
   
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
   * istream extractor for a UnitCell<3>::LatticeSystem.
   *
   * \param  in       input stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, 
                              UnitCell<3>::LatticeSystem& lattice);

   /**
   * ostream inserter for an UnitCell<3>::LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, 
                              UnitCell<3>::LatticeSystem lattice);

} 
#endif 
