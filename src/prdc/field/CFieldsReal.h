#ifndef PRDC_C_FIELDS_REAL_H
#define PRDC_C_FIELDS_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>   // member template
#include <pscf/math/IntVec.h>         // template with defaults

// Forward declarations
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;

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
   *  - A DArray of RFT (RField) containers holds valus of each field 
   *    on the nodes of a regular grid. This is accessed by the rgrid()
   *    and rgrid(int) member functions.
   *
   *  - A DArray of DArray<double> containers holds components of each
   *    field in a symmetry-adapted Fourier expansion (i.e., in basis
   *    format). This is accessed by the basis() and basis(int) member
   *    functions.
   *
   * The CFields container provides public non-const access to both field
   * representations, and does not automatically update one of these 
   * field representations when the other is modified. Maintenance of the 
   * intended relationship between the two data representations is instead
   * left as the responsibility of an object that owns this container. 
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
      * \param fieldIo  associated FIT (FieldIo) object
      */
      void setFieldIo(FIT const & fieldIo);

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
      * Set stored value of nMonomer.
      * 
      * This function may only be called once. The value of nMonomer must
      * be positive.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Allocate memory for fields in rgrid format.
      *
      * This function may only be called once.
      *
      * \param dimensions  dimensions of spatial mesh
      */
      void allocateRGrid(IntVec<D> const & dimensions);

      /**
      * Allocate or re-allocate memory for fields in basis format.
      *
      * This function may only be called once.
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
      /// \name Field Accessors (return by reference)
      ///@{

      /**
      * Get array of all fields in basis format (non-const).
      */
      DArray< DArray<double> > & basis();

      /**
      * Get array of all fields in basis format (const)
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> > const & basis() const;

      /**
      * Get the field for one monomer type in basis format (non-const).
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> & basis(int monomerId);

      /**
      * Get the field for one monomer type in basis format (const)
      *
      * \param monomerId integer monomer type index (0, ... ,nMonomer-1)
      */
      DArray<double> const & basis(int monomerId) const;

      /**
      * Get array of all fields in r-grid format (non-const).
      */
      DArray<RFT> & rgrid();

      /**
      * Get array of all fields in r-grid format (const).
      */
      DArray<RFT> const & rgrid() const;

      /**
      * Get field for one monomer type in r-grid format (non-const)
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RFT & rgrid(int monomerId);

      /**
      * Get field for one monomer type in r-grid format (const).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RFT const & rgrid(int monomerId) const;

      ///@}
      /// \name Field Output
      ///@{

      /**
      * Write fields to an input stream in symmetrized basis format.
      *
      * \param out  output stream to which to write fields
      */
      void writeBasis(std::ostream& out) const;

      /**
      * Write fields to a named file, in symmetrized basis format.
      *
      * \param filename  name of file to which to write fields
      */
      void writeBasis(std::string filename) const;

      /**
      * Writes fields to an input stream in real-space (r-grid) format.
      *
      * \param out  output stream to which to write fields
      */
      void writeRGrid(std::ostream& out) const;

      /**
      * Writes fields to a named file in real-space (r-grid) format.
      *
      * \param filename  name of file to which to write fields
      */
      void writeRGrid(std::string filename) const;

      ///@}
      /// \name Boolean Variable Queries
      ///@{

      /**
      * Has memory been allocated for fields in r-grid format?
      */
      bool isAllocatedRGrid() const;

      /**
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis() const;

      /**
      * Does this container have up-to-date fields?
      */
      bool hasData() const;

      /**
      * Are the fields invariant under elements of the space group?
      */
      bool isSymmetric() const;

      ///@}
      /// \name Boolean Variable Setters
      ///@{

      /**
      * Set the hasData flag.
      *
      * This should be set true when fields are set to those computed
      * from the current w fields, and false when any input to that
      * calculation changes.
      */
      void setHasData(bool hasData);

      /**
      * Set the isSymmetric flag.
      *
      * This should be set true if and only if the fields are known to
      * have been computed from symmetric w fields, and the basis
      * representation exists.
      */
      void setIsSymmetric(bool isSymmetric);

      ///@}

   protected:

      /**
      * Get associated FieldIo object (const reference).
      */
      FIT const & fieldIo() const;

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

      /*
      * Does this container hold up-to-date field data?
      */
      bool hasData_;

      /*
      * Are the fields symmetric?
      */
      bool isSymmetric_;

   };

   // Public inline member functions

   // Get array of all fields in basis format (non-const)
   template <int D, class RFT, class FIT> inline
   DArray< DArray<double> >& CFieldsReal<D,RFT,FIT>::basis()
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_;
   }

   // Get array of all fields in basis format (const)
   template <int D, class RFT, class FIT> inline
   DArray< DArray<double> > const & CFieldsReal<D,RFT,FIT>::basis() const
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_;
   }

   // Get one field in basis format (non-const)
   template <int D, class RFT, class FIT> inline
   DArray<double> & CFieldsReal<D,RFT,FIT>::basis(int id)
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_[id];
   }

   // Get one field in basis format (const)
   template <int D, class RFT, class FIT> inline
   DArray<double> const & CFieldsReal<D,RFT,FIT>::basis(int id)
   const
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_[id];
   }

   // Get all fields in r-grid format (non-const)
   template <int D, class RFT, class FIT> inline
   DArray<RFT>& CFieldsReal<D,RFT,FIT>::rgrid()
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_;
   }

   // Get all fields in r-grid format (const)
   template <int D, class RFT, class FIT> inline
   DArray<RFT> const & CFieldsReal<D,RFT,FIT>::rgrid() const
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_;
   }

   // Get one field in r-grid format (non-const)
   template <int D, class RFT, class FIT> inline
   RFT& CFieldsReal<D,RFT,FIT>::rgrid(int id)
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_[id];
   }

   // Get one field in r-grid format (const)
   template <int D, class RFT, class FIT> inline
   RFT const & CFieldsReal<D,RFT,FIT>::rgrid(int id) const
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_[id];
   }

   // Has memory been allocated for fields in r-grid format?
   template <int D, class RFT, class FIT> inline 
   bool CFieldsReal<D,RFT,FIT>::isAllocatedRGrid() const
   {  return isAllocatedRGrid_; }

   // Has memory been allocated for fields in basis format?
   template <int D, class RFT, class FIT> inline 
   bool CFieldsReal<D,RFT,FIT>::isAllocatedBasis() const
   {  return isAllocatedBasis_; }

   // Are the fields up-to-date?
   template <int D, class RFT, class FIT> inline 
   bool CFieldsReal<D,RFT,FIT>::hasData() const
   {  return hasData_; }

   // Are the fields symmetric under elements of the space group?
   template <int D, class RFT, class FIT> inline 
   bool CFieldsReal<D,RFT,FIT>::isSymmetric() const
   {  return isSymmetric_; }

   // Protected inline member function

   // Associated FieldIo object (const reference).
   template <int D, class RFT, class FIT>
   inline 
   FIT const & CFieldsReal<D,RFT,FIT>::fieldIo() const
   {
      UTIL_CHECK(fieldIoPtr_);
      return *fieldIoPtr_;
   }

} // namespace Prdc
} // namespace Pscf
#endif
