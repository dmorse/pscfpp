#ifndef PRDC_C_FIELDS_REAL_H
#define PRDC_C_FIELDS_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>      // member template

// Forward declarations
namespace Util {
  template <int D> class IntVec;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * A list of c fields stored in both basis and r-grid format.
   *
   * <b> Template parameters </b>: The template parameters represent:
   * 
   *     - D   : integer dimensionality of space, D=1,2, or 3
   *     - RFT : field type for r-grid data (e.g., RField<D>)
   *     - FIT : FieldIo type for field io operations (e.g., FieldIo<D>)
   * 
   * <b> Field Representations </b>: A CFieldsReal contains a list of
   * nMonomer fields that are each associated with a monomer type. The
   * fields may be stored in two different formats:
   *
   *  - A DArray of RFT containers holds valus of each field on
   *    the nodes of a regular grid. This is accessed by the rgrid()
   *    and rgrid(int) member functions.
   *
   *  - A DArray of DArray<double> containers holds components of each
   *    field in a symmetry-adapted Fourier expansion (i.e., in basis
   *    format). This is accessed by the basis() and basis(int) member
   *    functions.
   *
   * <b> Subclasses </b>: Partial specializations of CFieldsReal are
   * used as base classes for classes Prdc::CFieldsReal \<D \> and 
   * Rpg::CFieldsReal \<D\>:
   *
   *  - Subclass Prdc::CFieldsReal \<D\> is derived from a partial
   *    specialization of CFieldsReal with template parameters 
   *    RFT = Cpu::RField\<D\> and FIT = Prdc::FieldIo\<D\> , and is 
   *    used in the pscf_pc CPU program.
   *
   *  - Subclass Rpg::CFieldsReal \<D\> is derived from a partial
   *    specialization of CFieldsReal with template parameters 
   *    RFT = Cuda::RField \<D\> and FIT = Rpg::FieldIo \<D\> , and 
   *    is used in the pscf_pg GPU accelerated program.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class RFT, class FIT>
   class CFieldsReal 
   {

   public:

      /**
      * Constructor.
      */
      CFieldsReal();

      /**
      * Destructor.
      */
      ~CFieldsReal();

      /// \name Initialization and Memory Management
      ///@{

      /**
      * Create association with FIT (store pointer).
      *
      * \param fieldIo  associated FIT object
      */
      void setFieldIo(FIT const & fieldIo);

      /**
      * Set stored value of nMonomer.
      * 
      * May only be called once.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Set unit cell used when writing field files.
      *
      * This function creates a stored pointer to a UnitCell<D> that is
      * is used by the writeBasis and writeRGrid functions, which each
      * write the unit cell parameters from in this object to a field 
      * file header. This function may only be called once.
      *
      * \param cell  unit cell that is used by writeBasis and writeRGrid.
      */
      void setWriteUnitCell(UnitCell<D> const & cell);

      /**
      * Allocate or re-allocate memory for fields in rgrid format.
      *
      * \param dimensions  dimensions of spatial mesh
      */
      void allocateRGrid(IntVec<D> const & dimensions);

      /**
      * Allocate or re-allocate memory for fields in basis format.
      *
      * \param nBasis  number of basis functions 
      */
      void allocateBasis(int nBasis);

      /**
      * Allocate memory for both r-grid and basis field formats.
      *
      * This function may only be called once.
      *
      * \param nMonomer  number of monomer types
      * \param nBasis  number of basis functions 
      * \param dimensions  dimensions of spatial mesh
      */
      void allocate(int nMonomer, int nBasis, IntVec<D> const & dimensions);

      ///@}
      /// \name Field Mutators and Accessors (return by reference)
      ///@{

      /**
      * Get array of all fields in basis format (non-const).
      */
      DArray< DArray<double> > & basis();
      {  return basis_; }

      /**
      * Get array of all fields in basis format (const)
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> > const & basis() const
      {  return basis_; }

      /**
      * Get the field for one monomer type in basis format (non-const).
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> & basis(int monomerId)
      {  return basis_[monomerId]; }

      /**
      * Get the field for one monomer type in basis format (const)
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> const & basis(int monomerId) const
      {  return basis_[monomerId]; }

      /**
      * Get array of all fields in r-grid format (non-const).
      */
      DArray<RFT> & rgrid()
      {  return rgrid_; }

      /**
      * Get array of all fields in r-grid format (const).
      */
      DArray<RFT> const & rgrid() const
      {  return rgrid_; }

      /**
      * Get field for one monomer type in r-grid format (non-const)
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RFT & rgrid(int monomerId)
      {  return rgrid_[monomerId]; }

      /**
      * Get field for one monomer type in r-grid format (const).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RFT const & rgrid(int monomerId) const
      {  return rgrid_[monomerId]; }

      ///@}
      /// \name Boolean Queries
      ///@{

      /**
      * Has memory been allocated for fields in r-grid format?
      */
      bool isAllocatedRGrid() const
      {  return isAllocatedRGrid_; }

      /**
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis() const
      {  return isAllocatedBasis_; }

      ///@}

   protected:

      /**
      * Get associated FIT object (const reference).
      */
      FIT const & fieldIo() const;
      {
         UTIL_CHECK(fieldIoPtr_);
         return *fieldIoPtr_;
      }

   private:

      /*
      * Array of fields in symmetry-adapted basis format
      *
      * Element basis_[i] is an array that contains the components
      * of the field associated with monomer i, in a symmetry-adapted
      * Fourier basis expansion. 
      */
      DArray< DArray<double> > basis_;

      /*
      * Array of fields in real-space grid (r-grid) format
      *
      * Element basis_[i] is an RFT that contains values of the 
      * field associated with monomer i on the nodes of a regular mesh.
      */
      DArray<RFT> rgrid_;

      /*
      * Number of monomer types.
      */
      int nMonomer_;

      /*
      * Pointer to associated UnitCell<D> object.
      */
      UnitCell<D> const * writeUnitCellPtr_;

      /*
      * Pointer to associated FIT (FieldIo) object
      */
      FIT const * fieldIoPtr_;

      /*
      * Has memory been allocated for fields in r-grid format?
      */
      bool isAllocatedRGrid_;

      /*
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis_;

   };

} // namespace Prdc
} // namespace Pscf
#endif
