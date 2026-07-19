#ifndef RP_C_FIELDS_H
#define RP_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
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
   namespace Rp {
      template <int D, class T> class FieldIo;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * A list of c fields stored in both basis and r-grid format.
   *
   * <b> Template parameters </b>: The template parameters represent:
   *
   *   - D   : integer dimensionality of space (D=1, 2, or 3)
   *   - T   : a "Types" class (Cpp<D> or Rpg::Types<D>)
   *
   * <b> Field Representations </b>: A CFields container has a list of
   * nMonomer fields that are each associated with a monomer type. The
   * fields may be stored in two different formats:
   *
   *  - A DArray of RField containers holds valus of each field on
   *    the nodes of a regular grid. This is accessed by the rgrid()
   *    and rgrid(int) member functions.
   *
   *  - A DArray of DArray<double> containers holds components of each
   *    field in a symmetry-adapted Fourier expansion (i.e., in basis
   *    format). This is accessed by the basis() and basis(int) member
   *    functions.
   *
   * The CFields container provides public non-const (read/write) access 
   * to both field representations, and does not automatically update one 
   * of these field representations when the other is modified. 
   *
   * <b> Subclasses </b>: Specializations of this class template are used
   * as base classes for two closely analogous class templates, also named 
   * CFields, that are defined in Rpc and Rpg namespaces for use in the 
   * pscf_rpc and pscf_rpg programs, respectively.
   *
   * \ingroup Rp_Field_Module
   */
   template <int D, class T>
   class CFields
   {

   public:

      /**
      * Constructor.
      */
      CFields();

      /**
      * Destructor.
      */
      ~CFields() = default;

      /// \name Initialization and Memory Management
      ///@{

      /**
      * Create association with a FieldIo (store pointer).
      *
      * \param fieldIo  associated FieldIo object
      */
      void setFieldIo(FieldIo<D,T> const & fieldIo);

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
      DArray<typename T::RField> & rgrid();

      /**
      * Get array of all fields in r-grid format (const).
      */
      DArray<typename T::RField> const & rgrid() const;

      /**
      * Get field for one monomer type in r-grid format (non-const)
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      typename T::RField & rgrid(int monomerId);

      /**
      * Get field for one monomer type in r-grid format (const).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      typename T::RField const & rgrid(int monomerId) const;

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
      FieldIo<D,T> const & fieldIo() const;

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
      * Element basis_[i] is an typename T::RField that contains values of the
      * field associated with monomer i on the nodes of a regular mesh.
      */
      DArray<typename T::RField> rgrid_;

      /*
      * Number of monomer types.
      */
      int nMonomer_;

      /*
      * Pointer to associated UnitCell<D> object.
      */
      UnitCell<D> const * writeUnitCellPtr_;

      /*
      * Pointer to associated FieldIo<D,T> (FieldIo) object
      */
      FieldIo<D,T> const * fieldIoPtr_;

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
   template <int D, class T> inline
   DArray< DArray<double> >& CFields<D,T>::basis()
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_;
   }

   // Get array of all fields in basis format (const)
   template <int D, class T> inline
   DArray< DArray<double> > const & CFields<D,T>::basis() const
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_;
   }

   // Get one field in basis format (non-const)
   template <int D, class T> inline
   DArray<double> & CFields<D,T>::basis(int id)
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_[id];
   }

   // Get one field in basis format (const)
   template <int D, class T> inline
   DArray<double> const & CFields<D,T>::basis(int id)
   const
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_[id];
   }

   // Get all fields in r-grid format (non-const)
   template <int D, class T> inline
   DArray<typename T::RField>& CFields<D,T>::rgrid()
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_;
   }

   // Get all fields in r-grid format (const)
   template <int D, class T> inline
   DArray<typename T::RField> const & CFields<D,T>::rgrid() const
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_;
   }

   // Get one field in r-grid format (non-const)
   template <int D, class T> inline
   typename T::RField& CFields<D,T>::rgrid(int id)
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_[id];
   }

   // Get one field in r-grid format (const)
   template <int D, class T> inline
   typename T::RField const & CFields<D,T>::rgrid(int id) const
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_[id];
   }

   // Has memory been allocated for fields in r-grid format?
   template <int D, class T> inline
   bool CFields<D,T>::isAllocatedRGrid() const
   {  return isAllocatedRGrid_; }

   // Has memory been allocated for fields in basis format?
   template <int D, class T> inline
   bool CFields<D,T>::isAllocatedBasis() const
   {  return isAllocatedBasis_; }

   // Are the fields up-to-date?
   template <int D, class T> inline
   bool CFields<D,T>::hasData() const
   {  return hasData_; }

   // Are the fields symmetric under elements of the space group?
   template <int D, class T> inline
   bool CFields<D,T>::isSymmetric() const
   {  return isSymmetric_; }

   // Protected inline member function

   // Associated FieldIo object (const reference).
   template <int D, class T> inline
   FieldIo<D,T> const & CFields<D,T>::fieldIo() const
   {
      UTIL_CHECK(fieldIoPtr_);
      return *fieldIoPtr_;
   }

} // namespace Rp
} // namespace Pscf
#endif
