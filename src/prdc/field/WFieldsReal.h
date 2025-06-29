#ifndef PRDC_W_FIELDS_REAL_H
#define PRDC_W_FIELDS_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>              // member
#include <util/containers/DArray.h>        // member template

// Forward declarations
namespace Util {
   template <typename T> class Signal;
   template <> class Signal<void>;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {


   using namespace Util;

   /**
   * A container of input fields stored in both basis and r-grid format.
   * 
   * <b> Template parameters </b>: The template parameters represent:
   * 
   *     - D   : integer dimensionality of space, D=1,2, or 3
   *     - RFT : field type for r-grid data (e.g., RField<D>)
   *     - FIT : FieldIo type for field io operations (e.g., FieldIo<D>)
   * 
   * <b> Field Representations </b>: A WFieldsReal contains a list of
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
   * A WFieldsReal is designed to automatically update one of these
   * representations when the other is modified, when appropriate. A 
   * pointer to an associated FIT object is used for these conversions.
   *
   * The setBasis and readBasis functions allow the user to input new 
   * components in basis format, and both internally recompute the values 
   * in r-grid format.  The setRGrid and readRGrid functions allow the 
   * user to input the fields in r-grid format, and compute corresponding
   * components in basis format if and only if the user declares that the 
   * fields are known to be invariant under all symmetries of the space 
   * group. A boolean flag named isSymmetric is used to keep track of 
   * whether the current field is symmetric, and thus whether the basis 
   * format exists.
   *
   * <b> Subclasses </b>: Partial specializations of WFieldsReal are
   * used as base classes for classes Rpc::WFieldContainer \<D \> and 
   * Rpg::WFieldContainer \<D\>:
   *
   *  - Subclass Rpc::WFieldContainer \<D\> is derived from a partial
   *    specialization of WFieldsReal with template parameters 
   *    RFT = Cpu::RFT \<D\> and FIT = Rpc::FIT \<D\> , and is used in
   *    the pscf_pc CPU program.
   *
   *  - Subclass Rpg::WFieldContainer \<D\> is derived from a partial
   *    specialization of WFieldsReal with template parameters 
   *    RFT = Cuda::RFT \<D\> and FIT = Rpg::FIT \<D\> , and is used in
   *    the pscf_pg GPU accelerated program.
   *
   * <b> Signal </b>: A WFieldsReal owns an instance of class
   * Util::Signal<void> that notifies all observers whenever the fields
   * owned by the WFieldsReal are modified. This signal object may be 
   * accessed by reference using the signal() member function. The
   * Util::Signal<void>::addObserver function may used to add "observer"
   * objects and indicate a zero-parameter member function of each 
   * observer that will be called whenever the fields are modified.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class RFT, class FIT>
   class WFieldsReal 
   {

   public:

      /**
      * Constructor.
      */
      WFieldsReal();

      /**
      * Destructor.
      */
      ~WFieldsReal();

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
      * Set unit cell used when reading field files. 
      *
      * This function creates a stored pointer to a UnitCell<D> that is
      * is used by the readBasis and readRGrid functions, which reset the
      * unit cell parameters in this object to those read from the field 
      * file header. This function may only be called once.
      *
      * \param cell  unit cell that is modified by readBasis and readRGrid.
      */
      void setReadUnitCell(UnitCell<D>& cell);

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
      * Allocate memory for all fields.
      *
      * This function may only be called once.
      *
      * \param nMonomer  number of monomer types
      * \param nBasis  number of basis functions
      * \param dimensions  dimensions of spatial mesh
      */
      void allocate(int nMonomer, int nBasis, IntVec<D> const & dimensions);

      ///@}
      /// \name Field Modifiers
      ///@{

      /**
      * Set field component values, in symmetrized Fourier format.
      *
      * This function also computes and stores the corresponding r-grid
      * representation. On return, hasData and isSymmetric are both true.
      * 
      * The associated basis must be initialized on entry.  As needed, 
      * r-grid and/or basis fields may be allocated within this function.
      *
      * \param fields  array of new fields in basis format
      */
      void setBasis(DArray< DArray<double> > const & fields);

      /**
      * Set fields values in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that
      * the fields are known to be symmetric and so computes and stores
      * the corresponding basis components. If isSymmetric is false, it
      * only sets the values in the r-grid format.
      *
      * On return, hasData is true and the persistent isSymmetric flag
      * defined by the class is set to the value of the isSymmetric
      * input parameter.
      *
      * As needed, r-grid and/or basis fields may be allocated within this 
      * function. If the isSymmetric parameter is true, the a basis must 
      * be initialized prior to entry.
      *
      * \param fields  array of new fields in r-grid format
      * \param isSymmetric is this field symmetric under the space group?
      */
      void setRGrid(DArray< RFT > const & fields,
                    bool isSymmetric = false);

      /**
      * Read fields from an input stream in symmetrized basis format.
      *
      * This function also computes and stores the corresponding r-grid
      * representation. On return, hasData and isSymmetric are both true.
      *
      * As needed, r-grid and/or basis fields can be allocated within
      * this function, if not allocated on entry. An associated basis 
      * will be initialized if not initialized on entry.
      *
      * \param in  input stream from which to read fields
      */
      void readBasis(std::istream& in);

      /**
      * Read fields from a named file, in symmetrized basis format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      *
      * As needed, r-grid and/or basis fields may be allocated within
      * this function, if not allocated on entry. An associated basis
      * will be initialized if not initialized on entry.
      *
      * \param filename  file from which to read fields
      */
      void readBasis(std::string filename);

      /**
      * Reads fields from an input stream in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that
      * the fields are known to be symmetric and so computes and stores
      * the corresponding basis components. If isSymmetric is false, it
      * only sets the values in the r-grid format.
      *
      * On return, hasData is true and the persistent isSymmetric flag
      * defined by the class is set to the value of the isSymmetric
      * input parameter.
      *
      * As needed, r-grid and/or basis fields may be allocated within
      * this function, if not allocated on entry. An associated basis
      * will be initialized if not initialized on entry.
      *
      * \param in  input stream from which to read fields
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::istream& in, bool isSymmetric = false);

      /**
      * Reads fields from a named file in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that
      * the fields are known to be symmetric and so computes and stores
      * the corresponding basis components. If isSymmetric is false, it
      * only sets the values in the r-grid format.
      *
      * On return, hasData is true and the persistent isSymmetric flag
      * defined by the class is set to the value of the isSymmetric input
      * parameter.
      *
      * As needed, r-grid and/or basis fields may be allocated within
      * this function, if not allocated on entry. An associated basis
      * will be initialized if not initialized on entry.
      *
      * \param filename  file from which to read fields
      * \param isSymmetric  Is this field symmetric under the space group?
      */
      void readRGrid(std::string filename, bool isSymmetric = false);

      /**
      * Symmetrize r-grid fields, compute corresponding basis components.
      *
      * This function may be used after setting or reading w fields in
      * r-grid format that are known to be symmetric under the space
      * group to remove small deviations from symmetry and generate
      * basis components.
      *
      * The function symmetrizes the fields by converting from r-grid
      * to basis format and then back again, while also storing the
      * resulting basis components and setting isSymmetric() true.
      *
      * This function assumes that the current wFieldsRGrid fields
      * are known by the user to be symmetric, and does NOT check this.
      * Applying this function to fields that are not symmetric will
      * silently corrupt the fields.
      *
      * On entry, hasData() must be true and isSymmetric() must be false.
      * On exit, isSymmetric() is true.
      */
      void symmetrize();

      /**
      * Clear data stored in this object without deallocating.
      */
      void clear();

      /**
      * Get a signal that notifies observers of field modification.
      */
      Signal<void>& signal();

      ///@}
      /// \name Field Accessors (by const reference)
      ///@{

      /**
      * Get array of all fields in basis format.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> > const & basis() const;

      /**
      * Get the field for one monomer type in basis format.
      *
      * An Exception is thrown if isSymmetric is false.
      *
      * \param monomerId integer monomer type index (0,...,nMonomer-1)
      */
      DArray<double> const & basis(int monomerId) const;

      /**
      * Get array of all fields in r-space grid format.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< RFT > const & rgrid() const;

      /**
      * Get the field for one monomer type in r-space grid format.
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RFT const & rgrid(int monomerId) const;

      ///@}
      /// \name Boolean Queries
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
      * Has field data been set in either format?
      *
      * This flag is set true in setBasis and setRGrid.
      */
      bool hasData() const;

      /**
      * Are fields symmetric under all elements of the space group?
      *
      * A valid basis format exists if and only if isSymmetric is true.
      * This flag is set true when the fields are input in basis format
      * by the function setBasis or readBasis, or when they are set in 
      * grid format by calling the function setRGrid or readRGrid with 
      * function parameter isSymmetric == true.
      */
      bool isSymmetric() const;

      ///@}

   protected:

      /**
      * Get mesh dimensions in each direction, set on r-grid allocation.
      */
      IntVec<D> const & meshDimensions() const;

      /**
      * Get mesh size (number of grid points), set on r-grid allocation.
      */
      int meshSize() const;

      /**
      * Get number of basis functions, set on basis allocation.
      */
      int nBasis() const;

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      /**
      * Get associated FIT object (const reference).
      */
      FIT const & fieldIo() const;

   private:

      /*
      * Array of fields in symmetry-adapted basis format.
      *
      * Element basis_[i] is an array that contains the components
      * of the field associated with monomer i, in a symmetry-adapted
      * Fourier expansion.
      */
      DArray< DArray<double> > basis_;

      /*
      * Array of fields in real-space grid (r-grid) format.
      *
      * Element basis_[i] is an RFT that contains values of the
      * field associated with monomer i on the nodes of a regular mesh.
      */
      DArray< RFT > rgrid_;

      /*
      * Integer vector of grid dimensions.
      *
      * Element i is the number of grid points along direction i
      */
      IntVec<D> meshDimensions_;

      /*
      * Total number grid points (product of mesh dimensions).
      */
      int meshSize_;

      /*
      * Number of basis functions in symmetry-adapted basis.
      */
      int nBasis_;

      /*
      * Number of monomer types (number of fields).
      */
      int nMonomer_;

      /*
      * Pointer to unit cell modified by read functions.
      */
      UnitCell<D> * readUnitCellPtr_;

      /*
      * Pointer to unit cell access by write functions.
      */
      UnitCell<D> const * writeUnitCellPtr_;

      /*
      * Pointer to an associated FIT object.
      */
      FIT const * fieldIoPtr_;

      /*
      * Pointer to a Signal that is triggered by field modification.
      *
      * The Signal is constructed and owned by this field container.
      */
      Signal<void>* signalPtr_;
 
      /*
      * Has memory been allocated for fields in r-grid format?
      */
      bool isAllocatedRGrid_;

      /*
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis_;

      /*
      * Has field data been initialized ?
      */
      bool hasData_;

      /*
      * Are the fields symmetric under space group operations?
      *
      * Set true when fields are set using the symmetry adapted basis
      * format via function setBasis. False by otherwise.
      */
      bool isSymmetric_;

      /*
      *  Assign one RFT to another: lhs = rhs.
      *  
      *  \param lhs  left hand side of assignment
      *  \param rhs  right hand side of assignment
      */
      virtual void assignRField(RFT& lhs, RFT const & rhs) const;

   };

   // Inline member functions

   // Clear data stored in this object without deallocating
   template <int D, class RField, class FieldIo>
   inline void WFieldsReal<D,RField,FieldIo>::clear()
   {  hasData_ = false; }

   // Get array of all fields in basis format (const)
   template <int D, class RFT, class FIT>
   inline
   DArray< DArray<double> > const & 
   WFieldsReal<D,RFT,FIT>::basis() const
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_;
   }

   // Get one field in basis format (const)
   template <int D, class RFT, class FIT>
   inline
   DArray<double> const & WFieldsReal<D,RFT,FIT>::basis(int id)
   const
   {
      UTIL_ASSERT(isAllocatedBasis_);
      return basis_[id];
   }

   // Get all fields in r-grid format (const)
   template <int D, class RFT, class FIT>
   inline
   DArray< RFT > const &
   WFieldsReal<D,RFT,FIT>::rgrid() const
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_;
   }

   // Get one field in r-grid format (const)
   template <int D, class RFT, class FIT>
   inline
   RFT const & WFieldsReal<D,RFT,FIT>::rgrid(int id) const
   {
      UTIL_ASSERT(isAllocatedRGrid_);
      return rgrid_[id];
   }

   // Has memory been allocated for fields in r-grid format?
   template <int D, class RFT, class FIT>
   inline 
   bool WFieldsReal<D,RFT,FIT>::isAllocatedRGrid() const
   {  return isAllocatedRGrid_; }

   // Has memory been allocated for fields in basis format?
   template <int D, class RFT, class FIT>
   inline 
   bool WFieldsReal<D,RFT,FIT>::isAllocatedBasis() const
   {  return isAllocatedBasis_; }

   // Has field data been initialized ?
   template <int D, class RFT, class FIT>
   inline 
   bool WFieldsReal<D,RFT,FIT>::hasData() const
   {  return hasData_; }

   // Are the fields symmetric under space group operations?
   template <int D, class RFT, class FIT>
   inline 
   bool WFieldsReal<D,RFT,FIT>::isSymmetric() const
   {  return isSymmetric_; }

   // Protected inline member functions
   
   // Get mesh dimensions in each direction, set on r-grid allocation.
   template <int D, class RFT, class FIT>
   inline 
   IntVec<D> const & 
   WFieldsReal<D,RFT,FIT>::meshDimensions() const
   {  return meshDimensions_; }

   // Get mesh size (number of grid points), set on r-grid allocation.
   template <int D, class RFT, class FIT>
   inline 
   int WFieldsReal<D,RFT,FIT>::meshSize() const
   {  return meshSize_; }

   // Get number of basis functions, set on basis allocation.
   template <int D, class RFT, class FIT>
   inline 
   int WFieldsReal<D,RFT,FIT>::nBasis() const
   {  return nBasis_; }

   // Get number of monomer types.
   template <int D, class RFT, class FIT>
   inline 
   int WFieldsReal<D,RFT,FIT>::nMonomer() const
   {  return nMonomer_; }

   // Associated FieldIo object (const reference).
   template <int D, class RFT, class FIT>
   inline 
   FIT const & WFieldsReal<D,RFT,FIT>::fieldIo() const
   {
      UTIL_CHECK(fieldIoPtr_);
      return *fieldIoPtr_;
   }

} // namespace Prdc
} // namespace Pscf
#endif
