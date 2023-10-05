#ifndef PSPC_MASK_H
#define PSPC_MASK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>     // base class

#include <pspc/field/RField.h>             // member template parameter
#include <prdc/crystal/UnitCell.h>         // function parameter
#include <pscf/math/IntVec.h>              // function parameter
#include <util/containers/DArray.h>        // member template

namespace Pscf {
namespace Pspc {

   template <int D> class FieldIo;

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * A field to which the total density is constrained.
   *
   * A Mask<D> contains representations of a single field in two formats:
   * 
   *  - A DArray<double> that contains components of the field in a
   *    symmetry-adapted Fourier expansion (i.e., in basis format). This 
   *    is accessed by the basis() and basis(int) member functions.
   *
   *  - An RField<D> object contains values of the field on the nodes of a
   *    regular mesh. This is accessed by the rgrid() and rgrid(int) member 
   *    functions.
   *
   * A Mask is designed to automatically update one of these 
   * representations when the other is modified, when appropriate. 
   * A pointer to an associated FieldIo<D> is used for these conversions.
   * The setBasis function allows the user to input new components in
   * basis format and internally recomputes the values in r-grid format.
   * The setRgrid function allows the user to reset the field in 
   * r-grid format, but recomputes the components in basis format if 
   * and only if the user declares that the field are known to be
   * invariant under all symmetries of the space group. A boolean
   * flag named isSymmetric is used to keep track of whether the 
   * current field is symmetric, and thus whether the basis format 
   * exists.
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class Mask : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Mask();

      /**
      * Destructor.
      */
      ~Mask();

      /**
      * Create association with FieldIo (store pointer).
      */
      void setFieldIo(FieldIo<D> const & fieldIo);

      /**
      * Allocate memory for the field.
      *
      * A Mask<D> may only be allocated once. An Exception will
      * be thrown if this function is called more than once.
      *
      * \param nBasis  number of basis functions 
      * \param dimensions  dimensions of spatial mesh
      */
      void allocate(int nBasis, IntVec<D> const & dimensions);

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
      void setRGrid(RField<D> const & field, 
                    bool isSymmetric = false);

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

      /**
      * Get the field in basis format.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<double> const & basis() const;

      /**
      * Get array of all field in r-space grid format.
      *
      * The array capacity is equal to the number of monomer types.
      */
      RField<D> const & rgrid() const;

      /**
      * Volume fraction of the unit cell occupied by the polymers/solvents.
      * 
      * This value is equivalent to the spatial average of the mask, which is 
      * the q=0 coefficient of the discrete Fourier transform.
      * 
      * If hasData_ = true and isSymmetric_ = false (data only exists in 
      * rgrid format), then this object must be associated with a FieldIo
      * object in order to call phiTot(). In other cases, the FieldIo 
      * association is not necessary.
      */
      double phiTot() const;

      /**
      * Has memory been allocated?
      */
      bool isAllocated() const;

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

   private:

      /*
      * Components of field in symmetry-adapted basis format
      */
      DArray<double> basis_;

      /*
      * Field in real-space grid (r-grid) format
      */
      RField<D> rgrid_;

      /*
      * Pointer to associated FieldIo object
      */
      FieldIo<D> const * fieldIoPtr_;

      /*
      * Integer vector of grid dimensions.
      *
      * Element i is the number of grid points along direction i
      */
      IntVec<D> meshDimensions_;

      /*
      * Total number grid points (product of mesh dimensions)
      */
      int meshSize_;

      /*
      * Number of basis functions in symmetry-adapted basis
      */
      int nBasis_;

      /*
      * Number of monomer types (number of field)
      */
      int nMonomer_;

      /*
      * Has memory been allocated for field?
      */
      bool isAllocated_;

      /*
      * Has field data been initialized ?
      */
      bool hasData_;

      /*
      * Are the field symmetric under space group operations?
      *
      * Set true when field are set using the symmetry adapted basis
      * format via function setBasis. False by otherwise.
      */
      bool isSymmetric_;

   };

   // Inline member functions

   // Get array of all field in basis format (const)
   template <int D>
   inline
   DArray<double> const & Mask<D>::basis() const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_;
   }

   // Get all field in r-grid format (const)
   template <int D>
   inline
   RField<D> const & Mask<D>::rgrid() const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_;
   }

   // Has memory been allocated?
   template <int D>
   inline bool Mask<D>::isAllocated() const
   {  return isAllocated_; }

   // Have the field data been set?
   template <int D>
   inline bool Mask<D>::hasData() const
   {  return hasData_; }

   // Is the field symmetric under space group operations?
   template <int D>
   inline bool Mask<D>::isSymmetric() const
   {  return isSymmetric_; }

   #ifndef PSPC_W_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class Mask<1>;
   extern template class Mask<2>;
   extern template class Mask<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
