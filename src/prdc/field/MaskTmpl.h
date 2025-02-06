#ifndef PRDC_MASK_TMPL_H
#define PRDC_MASK_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>         // function parameter
#include <pscf/math/IntVec.h>              // function parameter
#include <util/containers/DArray.h>        // member template
#include <util/param/ParamComposite.h>     // base class

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Container for a field to which the total density is constrained.
   * 
   * MaskTmpl<D, FieldIo, RField> is a template class for Mask<D>. A
   * system that contains a Mask must satisfy a modified version of the
   * incompressibility constraint, in which the sum of the concentration
   * fields of all species must be equal to the Mask field. The Mask takes
   * values from 0 to 1 everywhere. A system without a Mask is equivalent
   * to a system in which the mask is equal to 1 at all points.
   * 
   * A Mask<D> contains representations of this field in two formats:
   * 
   *  - A DArray<double> that contains components of the field in a
   *    symmetry-adapted Fourier expansion (i.e., in basis format). This 
   *    is accessed by the basis() member function.
   *
   *  - An RField object (where RField is a template class) contains values 
   *    of the field on the nodes of a regular mesh. This is accessed by 
   *    the rgrid() member function.
   *
   * A Mask is designed to automatically update one of these 
   * representations when the other is modified, when appropriate. 
   * A pointer to an associated FieldIo (another template class) is used 
   * for these conversions. The FieldIo class that is used to instantiate
   * this template should be a subclass of Prdc::FieldIoReal, and the 
   * RField class used to instantiate this template should be compatible
   * with the FieldIo.
   * 
   * The setBasis function allows the user to input new components in
   * basis format and internally recomputes the values in r-grid format.
   * The setRgrid function allows the user to reset the field in r-grid
   * format, but recomputes the components in basis format if and only
   * if the user explicitly declares that the field are known to be
   * invariant under all symmetries of the space group. A boolean flag
   * named isSymmetric is used to keep track of whether the current
   * field is symmetric, and thus whether the basis format exists.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, typename FieldIo, typename RField>
   class MaskTmpl : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      MaskTmpl();

      /**
      * Destructor.
      */
      ~MaskTmpl();

      /// \name Initialization and Memory Management
      ///@{

      /**
      * Create association with FieldIo (store pointer).
      */
      void setFieldIo(FieldIo const & fieldIo);

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
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
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
      * \param field  new field in r-grid format
      * \param isSymmetric is this field symmetric under the space group?
      */
      void setRGrid(RField const & field, bool isSymmetric = false);

      /**
      * Read field from input stream in symmetrized Fourier format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      * 
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      *
      * \param in  input stream from which to read field
      * \param unitCell  associated crystallographic unit cell
      */
      void readBasis(std::istream& in, UnitCell<D>& unitCell);

      /**
      * Read field from file in symmetrized Fourier format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      * 
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      *
      * \param filename  file from which to read field
      * \param unitCell  associated crystallographic unit cell
      */
      void readBasis(std::string filename, UnitCell<D>& unitCell);

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
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      * 
      * \param in  input stream from which to read field
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::istream& in, UnitCell<D>& unitCell,
                     bool isSymmetric = false);

      /**
      * Reads field from a file in real-space (r-grid) format.
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
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      * 
      * \param filename  file from which to read field
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::string filename, UnitCell<D>& unitCell,
                     bool isSymmetric = false);

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
      RField const & rgrid() const;

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
      * Calculate the average value of the rgrid_ member.
      * 
      * Must be implemented by subclasses, because this calculation
      * depends on the specific implementation of the RField type.
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
      RField rgrid_;

      /**
      * Pointer to associated FieldIo object
      */
      FieldIo const * fieldIoPtr_;

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

      /**
      * Number of monomer types (number of field)
      */
      int nMonomer_;

      /**
      * Has memory been allocated for field in basis format?
      */
      bool isAllocatedBasis_;

      /**
      * Has memory been allocated for field in rgrid format?
      */
      bool isAllocatedRGrid_;

      /**
      * Has field data been initialized ?
      */
      bool hasData_;

      /**
      * Are the field symmetric under space group operations?
      *
      * Set true when field are set using the symmetry adapted basis
      * format via function setBasis. False by otherwise.
      */
      bool isSymmetric_;

   };

   // Inline member functions

   // Get field in basis format (const)
   template <int D, typename FieldIo, typename RField>
   inline DArray<double> const & MaskTmpl<D, FieldIo, RField>::basis() const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_;
   }

   // Get field in r-grid format (const)
   template <int D, typename FieldIo, typename RField>
   inline RField const & MaskTmpl<D, FieldIo, RField>::rgrid() const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_;
   }

   // Has memory been allocated in basis format?
   template <int D, typename FieldIo, typename RField>
   inline bool MaskTmpl<D, FieldIo, RField>::isAllocatedBasis() const
   {  return isAllocatedBasis_; }

   // Has memory been allocated in rgrid format?
   template <int D, typename FieldIo, typename RField>
   inline bool MaskTmpl<D, FieldIo, RField>::isAllocatedRGrid() const
   {  return isAllocatedRGrid_; }

   // Have the field data been set?
   template <int D, typename FieldIo, typename RField>
   inline bool MaskTmpl<D, FieldIo, RField>::hasData() const
   {  return hasData_; }

   // Is the field symmetric under space group operations?
   template <int D, typename FieldIo, typename RField>
   inline bool MaskTmpl<D, FieldIo, RField>::isSymmetric() const
   {  return isSymmetric_; }

} // namespace Prdc
} // namespace Pscf
#endif
