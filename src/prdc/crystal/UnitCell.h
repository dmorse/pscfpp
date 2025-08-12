#ifndef PRDC_UNIT_CELL_H
#define PRDC_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCellBase.h"
#include <util/format/Dbl.h>
 
#include <iostream>
#include <iomanip>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base template for UnitCell<D> classes, D=1, 2 or 3.
   *
   * Explicit specializations are defined for D=1, 2, and 3. In each 
   * case, class UnitCell<D> is derived from UnitCellBase<D> and defines 
   * an enumeration UnitCell<D>::LatticeSystem of the possible types of 
   * Bravais lattice systems in D-dimensional space.  Iostream inserter 
   * (<<) and extractor (>>) operators are defined for each
   * UnitCell<D>::LatticeSystem enumeration. These operators define a 
   * standard text representation of each allowed unit cell type, as
   * described \ref scft_unitcell_page "here".
   *
   * Each explicit specialization UnitCell<D> has a member variable of 
   * type UnitCell<D>::LatticeSystem, the value of which is returned by
   * a function named lattice().  The value of the lattice variable is 
   * initialized to an enumeration value named Null, which denotes 
   * unknown or uninitialized. This variable may only be set once. Any
   * attempt to reset the lattice type to a different non-null after it 
   * has been set will cause an Exception to be thrown and program 
   * execution to halt. Each UnitCell<D> object may thus only be used 
   * to represent one lattice type during its lifetime. 
   *
   * A template defines iostream inserter (<<) and extractor (>>) 
   * operators for each explicit specializations of UnitCell<D>, allowing 
   * a UnitCell<D> to be read from or written to file like a primitive 
   * variable. The text representation for a UnitCell<D> contains a text
   * representation of the LatticeSystem<D> enumeration (i.e., the unit
   * cell type) and a list of one or more unit cell parameters (lengths
   * and angles), as described \ref scft_unitcell_page "here".
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   class UnitCell : public UnitCellBase<D>
   {};

   // Function template declarations 

   /**
   * istream input extractor for a UnitCell<D>.
   *
   * If cell.lattice() is not Null, the lattice type read from stream
   * must match the existing lattice type, or an Exception is thrown.
   *
   * \param  in  input stream
   * \param  cell  UnitCell<D> to be read
   * \return modified input stream
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   std::istream& operator >> (std::istream& in, UnitCell<D>& cell);

   /**
   * ostream output inserter for a UnitCell<D>.
   *
   * \param out  output stream
   * \param  cell  UnitCell<D> to be written
   * \return modified output stream
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, UnitCell<D> const& cell);

   /**
   * Serialize to/from an archive.
   *
   * \param ar  input or output archive
   * \param cell  UnitCell<D> object to be serialized
   * \param version  archive version id
   * \ingroup Prdc_Crystal_Module
   */
   template <class Archive, int D>
   void 
   serialize(Archive& ar, UnitCell<D>& cell, const unsigned int version);

   // 1D Unit Cell

   /**
   * 1D crystal unit cell.
   *
   * \ingroup Prdc_Crystal_Module
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

      /**
      *  Assignment operator. 
      * 
      * Copies lattice type and unit cell parameters, and recomputes
      * all data that depends on the unit cell parameters, but does
      * not copy any associated Signal<> object.
      *      
      * \param other UnitCell<1> object to be cloned.
      */
      UnitCell<1>& operator = (const UnitCell<1>& other);

      /**
      * Set the lattice system, but not unit cell parameters.
      *
      * Upon return, values of lattice and nParameter are set. 
      *
      * If cell.lattice() is not Null, the new lattice type must match
      * match the existing lattice type, or an Exception is thrown. 
      * The lattice type may thus only be set once.
      *
      * \param lattice  lattice system enumeration value
      */
      void set(UnitCell<1>::LatticeSystem lattice);

      /**
      * Set the unit cell state (lattice system and parameters).
      *
      * This function calls set(UnitCell<1>::Lattice System) internally.
      * If cell.lattice() is not Null, the new lattice type must match
      * match the existing lattice type, or an Exception is thrown. 
      *
      * \param lattice  lattice system enumeration value
      * \param parameters  array of unit cell parameters
      */
      void set(UnitCell<1>::LatticeSystem lattice, 
               FSArray<double, 6> const & parameters);

      /**
      * Return lattice system enumeration value.
      *
      * This value is initialized to Null during construction.
      */
      LatticeSystem lattice() const
      {  return lattice_; }

      /** 
      * Get the generalized volume (i.e., length) of the unit cell,
      */
      double volume() const;

      using UnitCellBase<1>::isInitialized;

   private:

      // Lattice type (lamellar or Null)
      LatticeSystem lattice_;

      // Set number of parameters required to describe current lattice type
      void setNParameter();

      // Set all internal data after setting parameter values
      void setBasis();

      // Private and unimplemented to prevent copy construction.
      UnitCell(UnitCell<1> const &);

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
   * \ingroup Prdc_Crystal_Module
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<1>::LatticeSystem& lattice);

   /**
   * ostream inserter for a 1D UnitCell<1>::LatticeSystem.
   *
   * \param out  output stream
   * \param lattice  UnitCell<1>::LatticeSystem to be written
   * \return modified output stream
   * \ingroup Prdc_Crystal_Module
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<1>::LatticeSystem lattice);

   /**
   * Serialize a UnitCell<1>::LatticeSystem enumeration value
   *
   * \param ar  archive
   * \param lattice  enumeration data to be serialized
   * \param version  version id
   */ 
   template <class Archive>
   inline 
   void serialize(Archive& ar, UnitCell<1>::LatticeSystem& lattice, 
                  const unsigned int version)
   {  serializeEnum(ar, lattice, version); }


   // 2D Unit Cell

   /**
   * 2D crystal unit cell.
   *
   * \ingroup Prdc_Crystal_Module
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

      /**
      *  Assignment operator. 
      *      
      * Copies lattice type and unit cell parameters, and recomputes
      * all data that depends on the unit cell parameters, but does
      * not copy any associated Signal<> object.
      *      
      * \param other UnitCell<2> object to be cloned.
      */
      UnitCell<2>& operator = (const UnitCell<2>& other);

      /**
      * Set the lattice system, but not unit cell parameters.
      *
      * Upon return, values of lattice and nParameter are set. 
      *
      * If cell.lattice() is not Null, the new lattice type must match
      * match the existing lattice type, or an Exception is thrown. 
      * The lattice type may thus only be set once.
      *
      * \param lattice  lattice system enumeration value
      */
      void set(UnitCell<2>::LatticeSystem lattice);

      /**
      * Set the unit cell state (lattice system and parameters).
      *
      * This function calls set(UnitCell<2>::LatticeSystem) internally.
      * If cell.lattice() is not Null, the new lattice type must match
      * match the existing lattice type, or an Exception is thrown. 
      * The lattice type may thus only be set once.
      *
      * \param lattice  lattice system enumeration value
      * \param parameters  array of unit cell parameters
      */
      void set(UnitCell<2>::LatticeSystem lattice, 
               FSArray<double, 6> const & parameters);

      /**
      * Return lattice system enumeration value.
      *
      * This value is initialized to Null during construction.
      */
      LatticeSystem lattice() const
      {  return lattice_; }

      /** 
      * Get the generalized volume (i.e., area) of the unit cell,
      */
      double volume() const;

   private:

      // Lattice system (square, rectangular, etc.)
      LatticeSystem lattice_;

      // Set number of parameters required to describe current lattice type
      void setNParameter();

      // Set all internal data after setting parameter values
      void setBasis();

      // Private and unimplemented to prevent copy construction.
      UnitCell(UnitCell<2> const &);

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
   * \ingroup Prdc_Crystal_Module
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
   * Serialize a UnitCell<2>::LatticeSystem enumeration value
   *
   * \param ar  archive
   * \param lattice  enumeration data to be serialized
   * \param version  version id
   */ 
   template <class Archive>
   inline 
   void serialize(Archive& ar, UnitCell<2>::LatticeSystem& lattice, 
                  const unsigned int version)
   {  serializeEnum(ar, lattice, version); }


   // 3D crystal unit cell

   /**
   * 3D crystal unit cell.
   *
   * \ingroup Prdc_Crystal_Module
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

      /**
      *  Assignment operator. 
      *      
      * Copies lattice type and unit cell parameters, and recomputes
      * all data that depends on the unit cell parameters, but does
      * not copy any associated Signal<> object.
      *      
      * \param other UnitCell<3> object to be cloned.
      */
      UnitCell<3>& operator = (const UnitCell<3>& other);

      /**
      * Set the lattice system, but not unit cell parameters.
      *
      * Upon return, values of lattice and nParameter are set. 
      *
      * If cell.lattice() is not Null, the new lattice type must match
      * match the existing lattice type, or an Exception is thrown. 
      * The lattice type may thus only be set once.
      *
      * \param lattice  lattice system enumeration value
      */
      void set(UnitCell<3>::LatticeSystem lattice);

      /**
      * Set the unit cell state (lattice type and cell parameters).
      *
      * This function calls set(UnitCell<3>::LatticeSystem) internally.
      * If cell.lattice() is not Null, the new lattice type must match
      * match the existing lattice type, or an Exception is thrown. 
      * The lattice type may thus only be set once.
      *
      * \param lattice  lattice system enumeration value
      * \param parameters  array of unit cell parameters
      */
      void set(UnitCell<3>::LatticeSystem lattice, 
               FSArray<double, 6> const & parameters);

      /**
      * Return lattice system enumeration value.
      *
      * This value is initialized to Null during construction.
      */
      LatticeSystem lattice() const
      {  return lattice_; }

      /** 
      * Get the volume of the unit cell,
      */
      double volume() const;

   private:

      LatticeSystem lattice_;

      // Set number of parameters required to describe current lattice type
      void setNParameter();

      // Set all internal data after setting parameter values
      void setBasis();

      // Private and unimplemented to prevent copy construction.
      UnitCell(UnitCell<3> const &);

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
   * \ingroup Prdc_Crystal_Module
   */
   std::istream& operator >> (std::istream& in,
                              UnitCell<3>::LatticeSystem& lattice);

   /**
   * ostream inserter for an 3D UnitCell<3>::LatticeSystem.
   *
   * \param  out      output stream
   * \param  lattice  UnitCell<3>::LatticeSystem to be written
   * \return modified output stream
   * \ingroup Prdc_Crystal_Module
   */
   std::ostream& operator << (std::ostream& out,
                              UnitCell<3>::LatticeSystem lattice);

   // General Function Templates

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

   /**
   * Serialize a UnitCell<3>::LatticeSystem enumeration value
   *
   * \param ar  archive
   * \param lattice  enumeration data to be serialized
   * \param version  version id
   */ 
   template <class Archive>
   inline 
   void serialize(Archive& ar, UnitCell<3>::LatticeSystem& lattice, 
                  const unsigned int version)
   {  serializeEnum(ar, lattice, version); }


   // UnitCell stream IO operator template definitions

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

   #ifndef PRDC_UNIT_CELL_TPP
   // Suppress implicit instantiation

   extern template class UnitCell<1>;
   extern template class UnitCell<2>;
   extern template class UnitCell<3>;

   extern std::ostream& operator << (std::ostream& out, 
                                     UnitCell<1>::LatticeSystem lattice);
   extern std::ostream& operator << (std::ostream& out, 
                                     UnitCell<2>::LatticeSystem lattice);
   extern std::ostream& operator << (std::ostream& out, 
                                     UnitCell<3>::LatticeSystem lattice);

   extern std::istream& operator >> (std::istream& in, 
                                UnitCell<1>::LatticeSystem const& lattice);
   extern std::istream& operator >> (std::istream& in, 
                                UnitCell<2>::LatticeSystem const& lattice);
   extern std::istream& operator >> (std::istream& in, 
                                UnitCell<3>::LatticeSystem const& lattice);
   #endif

}
}
#endif
