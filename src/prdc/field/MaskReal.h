#ifndef PRDC_MASK_TMPL_H
#define PRDC_MASK_TMPL_H

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
   * Container for a field to which the total density is constrained.
   * 
   * A system that contains a Mask must satisfy a modified version of the
   * incompressibility constraint, in which the sum of the concentration
   * (or volume fraction) fields of all monomer types must be equal to the 
   * Mask field. The Mask field takes values in the range [0, 1] everywhere. 
   * A system without a Mask is equivalent to one in which the mask field 
   * is equal to 1 at all points in the unit cell.
   * 
   * <b> Field representations </b>: A Mask \<D\> contains representations 
   * of the mask field in two formats:
   * 
   *  - An RFT object (where RFT is a template type parameter) contains
   *    values of the field on the nodes of a regular mesh. This is
   *    accessed by the rgrid() member function.
   *
   *  - A DArray \<double\> that contains components of the field in a
   *    symmetry-adapted Fourier expansion (i.e., in basis format). This 
   *    is accessed by the basis() member function.
   *
   * A Mask is designed to automatically update one of these
   * representations when the other is modified, when appropriate. 
   * A pointer to an associated FIT (another template parameter) is
   * used for these conversions. The FieldIo class that is used to 
   * instantiate this template should be a subclass of Prdc::FieldIoReal.
   * 
   * The setBasis and readBasis functions allow the user to input field
   * components in basis format, and both internally recompute the values 
   * in r-grid format.  The setRGrid and readRGrid functions allows the 
   * user to input the field in r-grid format, and both recompute the 
   * components in basis format if and only if the user explicitly declares 
   * that the field is known to be invariant under all symmetries of the 
   * space group. A boolean member variable named isSymmetric is used to 
   * keep track of whether the current mask field is symmetric, and thus 
   * whether the symmetry-adapted basis representation exists.
   *
   * <b> Subclasses </b>: Partial specializations of the template
   * MaskReal \<D, RFT, FIT\> are used as base classes for the class 
   * templates Rpc::Mask \<D \> and Rpg::Mask \<D\> that are used by
   * pscf_pc and pscf_pg, respectively.
   *
   * <b> Signal </b>: A MaskReal owns an instance of class
   * Util::Signal<void> that notifies all observers whenever the field
   * owned by the MaskReal is modified. This Signal object may be 
   * accessed by reference using the signal() member function. The
   * Util::Signal<void>::addObserver function may used to add observer 
   * objects and indicate a zero-parameter member function of each 
   * observer that will be called whenever the field is modified.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class RFT, class FIT>
   class MaskReal 
   {

   public:

      /**
      * Constructor.
      */
      MaskReal();

      /**
      * Destructor.
      */
      ~MaskReal();

      /// \name Initialization and Memory Management
      ///@{

      /**
      * Create association with FieldIo (store pointer).
      * 
      * \param fieldIo  associated FieldIo object
      */
      void setFieldIo(FIT const & fieldIo);

      /**
      * Set unit cell used when reading a mask field file.
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
      * Set unit cell used when writing a mask field file.
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
      * Allocate memory for the field in basis format.
      *
      * An Exception will be thrown if this is called more than once.
      *
      * \param nBasis  number of basis functions 
      */
      void allocateBasis(int nBasis);

      /**
      * Allocate memory for the field in rgrid format.
      *
      * An Exception will be thrown if this is called more than once.
      *
      * \param dimensions  dimensions of spatial mesh
      */
      void allocateRGrid(IntVec<D> const & dimensions);

      ///@}
      /// \name Field Mutators
      ///@{

      /**
      * Set field component values, in symmetrized Fourier format.
      *
      * This function also computes and stores the corresponding r-grid
      * representation. On return, hasData and isSymmetric are both true.
      *
      * The associated basis must be initialized on entry. As needed,
      * r-grid and/or basis fields may be allocated within this function.
      *
      * \param field  components of field in basis format
      */
      void setBasis(DArray<double> const & field);

      /**
      * Set field values in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that 
      * the field are known to be symmetric and so computes and stores
      * the corresponding basis components. If isSymmetric is false, it
      * only sets the values in the r-grid format.
      * 
      * On return, hasData is true and the persistent isSymmetric flag 
      * defined by the class is set to the value of the isSymmetric 
      * input parameter.
      *
      * As needed, r-grid and/or basis fields may be allocated within this
      * function. If the isSymmetric parameter is true, then a basis must
      * be initialized prior to entry.
      * 
      * \param field  new field in r-grid format
      * \param isSymmetric is this field symmetric under the space group?
      */
      void setRGrid(RFT const & field, bool isSymmetric = false);

      /**
      * Read field from input stream in symmetrized basis format.
      *
      * This function also computes and stores the corresponding r-grid
      * representation. On return, hasData and isSymmetric are both true.
      * 
      * This object must already be allocated and associated with a
      * a FieldIo object to run this function.
      *
      * As needed, r-grid and/or basis fields may be allocated within
      * this functions, if not allocated on entry. An associated basis
      * will be initialized if it is not initialized on entry.
      *
      * \param in  input stream from which to read field
      */
      void readBasis(std::istream& in);

      /**
      * Read field from a named file, in symmetrized basis format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      * 
      * As needed, r-grid and/or basis fields may be allocated within
      * this functions, if not allocated on entry. An associated basis
      * will be initialized if it is not initialized on entry.
      *
      * \param filename  file from which to read field
      */
      void readBasis(std::string filename);

      /**
      * Reads field from an input stream in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that 
      * the field is known to be symmetric and so computes and stores
      * the corresponding basis format. If isSymmetric is false, it
      * only sets the values in the r-grid format.
      * 
      * On return, hasData is true and the persistent isSymmetric flag 
      * defined by the class is set to the value of the isSymmetric 
      * input parameter.
      * 
      * As needed, r-grid and/or basis fields may be allocated within
      * this functions, if not allocated on entry. An associated basis
      * will be initialized if it is not initialized on entry.
      * 
      * \param in  input stream from which to read field
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::istream& in, bool isSymmetric = false);

      /**
      * Reads field from a named file, in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that 
      * the field is known to be symmetric and so computes and stores
      * the corresponding basis format. If isSymmetric is false, it
      * only sets the values in the r-grid format.
      * 
      * On return, hasData is true and the persistent isSymmetric flag 
      * defined by the class is set to the value of the isSymmetric 
      * input parameter.
      * 
      * As needed, r-grid and/or basis fields may be allocated within
      * this functions, if not allocated on entry. An associated basis
      * will be initialized if it is not initialized on entry.
      * 
      * \param filename  file from which to read field
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::string filename, bool isSymmetric = false);

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
      /// \name Field Accessors 
      ///@{

      /**
      * Get the field in basis format.
      */
      DArray<double> const & basis() const;

      /**
      * Get the field in r-grid format.
      */
      RFT const & rgrid() const;

      /**
      * Return the volume fraction of unit cell occupied by material.
      * 
      * This value is equivalent to the spatial average of the mask, 
      * which is the q=0 coefficient of the discrete Fourier transform.
      * 
      * If hasData == true and isSymmetric == false (i.e., data only 
      * exists in rgrid format), then this object must be associated 
      * with a FieldIo object in order to call phiTot(). In other 
      * cases, the FieldIo association is not necessary.
      */
      double phiTot() const;

      /**
      * Get a signal that notifies observers of field modification.
      *
      * The Signal<void>::notify method is called by within all member 
      * functions that modify the field, to notify observers of this 
      * event. The Signal<void>::addObserver function may be applied
      * to the Signal object returned by this function to add one or
      * more observers.
      */
      Signal<void>& signal();

      ///@}
      /// \name Boolean Queries
      ///@{

      /**
      * Has memory been allocated in basis format?
      */
      bool isAllocatedBasis() const;

      /**
      * Has memory been allocated in rgrid format?
      */
      bool isAllocatedRGrid() const;

      /**
      * Has field data been set in either format?
      *
      * This flag is set true in setBasis and setRGrid.
      */
      bool hasData() const;

      /**
      * Are field symmetric under all elements of the space group?
      *
      * A valid basis format exists if and only if isSymmetric is true.
      * This flat is set true if the field were input in basis format 
      * by the function setBasis, or if they were set in grid format
      * by the function setRGrid but isSymmetric was set true.
      */
      bool isSymmetric() const;

      ///@}

   protected:

      /**
      * Mesh dimensions in each direction, set by allocation.
      */
      IntVec<D> const & meshDimensions() const
      {  return meshDimensions_; }

      /**
      * Mesh size (number of grid points), set by allocation.
      */
      int meshSize() const
      {  return meshSize_; }

      /**
      * Number of basis functions, set by allocation.
      */
      int nBasis() const
      {  return nBasis_; }

      /**
      * Associated FieldIo object (const reference).
      */
      FIT const & fieldIo() const
      {
         UTIL_CHECK(fieldIoPtr_);
         return *fieldIoPtr_;
      }

      /**
      * Calculate the average value of the rgrid_ member.
      * 
      * Must be implemented by subclasses, because this calculation
      * depends on the specific implementation of the RFT type.
      */
      virtual double rGridAverage() const = 0;

   private:

      /**
      * Components of field in symmetry-adapted basis format
      */
      DArray<double> basis_;

      /**
      * Field in real-space grid (r-grid) format
      */
      RFT rgrid_;

      /**
      * Integer vector of grid dimensions.
      *
      * Element i is the number of grid points along direction i
      */
      IntVec<D> meshDimensions_;

      /**
      * Total number grid points (product of mesh dimensions)
      */
      int meshSize_;

      /**
      * Number of basis functions in symmetry-adapted basis
      */
      int nBasis_;

      /*
      * Pointer to unit cell modified by read functions.
      */
      UnitCell<D> * readUnitCellPtr_;

      /*
      * Pointer to unit cell access by write functions.
      */
      UnitCell<D> const * writeUnitCellPtr_;

      /**
      * Pointer to associated FieldIo object.
      */
      FIT const * fieldIoPtr_;

      /*
      * Pointer to a Signal that is triggered by field modification.
      */
      Signal<void>* signalPtr_;

      /**
      * Has memory been allocated for field in basis format?
      */
      bool isAllocatedBasis_;

      /**
      * Has memory been allocated for field in rgrid format?
      */
      bool isAllocatedRGrid_;

      /**
      * Has field data been initialized?
      */
      bool hasData_;

      /**
      * Is the mask field symmetric under space group operations?
      *
      * Set true when field are set using the symmetry adapted basis
      * format via function setBasis. False otherwise.
      */
      bool isSymmetric_;

   };

   // Inline member functions

   // Get field in basis format (const)
   template <int D, class RFT, class FIT>
   inline DArray<double> const & MaskReal<D,RFT,FIT>::basis() const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_;
   }

   // Get field in r-grid format (const)
   template <int D, class RFT, class FIT>
   inline RFT const & MaskReal<D,RFT,FIT>::rgrid() const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_;
   }

   // Has memory been allocated in basis format?
   template <int D, class RFT, class FIT>
   inline bool MaskReal<D,RFT,FIT>::isAllocatedBasis() const
   {  return isAllocatedBasis_; }

   // Has memory been allocated in rgrid format?
   template <int D, class RFT, class FIT>
   inline bool MaskReal<D,RFT,FIT>::isAllocatedRGrid() const
   {  return isAllocatedRGrid_; }

   // Have the field data been set?
   template <int D, class RFT, class FIT>
   inline bool MaskReal<D,RFT,FIT>::hasData() const
   {  return hasData_; }

   // Is the field symmetric under space group operations?
   template <int D, class RFT, class FIT>
   inline bool MaskReal<D,RFT,FIT>::isSymmetric() const
   {  return isSymmetric_; }

} // namespace Prdc
} // namespace Pscf
#endif
