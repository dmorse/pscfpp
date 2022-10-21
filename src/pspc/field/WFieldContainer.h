#ifndef PSPC_W_FIELD_CONTAINER_H
#define PSPC_W_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>     // base class

#include <pscf/math/IntVec.h>              // function parameter
#include <pscf/crystal/UnitCell.h>         // function parameter
#include <pspc/field/RField.h>             // member template parameter
#include <util/containers/DArray.h>        // member template

namespace Pscf {
namespace Pspc {

   template <int D> class FieldIo;

   using namespace Util;

   /**
   * A list of fields stored in both basis and r-grid format.
   *
   * A WFieldContainer<D> contains representations of a list of nMonomer
   * fields that are associated with different monomer types in two 
   * different related formats:
   * 
   *  - A DArray of DArray<double> containers holds components of each
   *    field in a symmetry-adapted Fourier expansion (i.e., in basis 
   *    format). This is accessed by the basis() and basis(int) 
   *    member functions.
   *
   *  - A DArray of RField<D> containers holds valus of each field on
   *    the nodes of a regular grid. This is accessed by the rgrid()
   *    and rgrid(int) member functions.
   *
   * A WFieldContainer is designed to automatically update one of these
   * representations when the other is modified, when appropriate. 
   * A pointer to an associated FieldIo<D> is used for these conversions.
   * The setBasis function allows the user to input new components in
   * basis format and internally recomputes the values in r-grid format.
   * The setRgrid function allows the user to reset the fields in 
   * r-grid format, but recomputes the components in basis format if 
   * and only if the user declares that the fields are known to be
   * invariant under all symmetries of the space group. A boolean
   * flag named isSymmetric is used to keep track of whether the 
   * current field is symmetric, and thus whether the basis format 
   * exists.
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class WFieldContainer : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      WFieldContainer();

      /**
      * Destructor.
      */
      ~WFieldContainer();

      /**
      * Create association with FieldIo (store pointer).
      */
      void setFieldIo(FieldIo<D> const & fieldIo);

      /**
      * Set stored value of nMonomer.
      * 
      * May only be called once.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Allocate or re-allocate memory for fields in rgrid format.
      *
      * \param dimensions  dimensions of spatial mesh
      */
      void allocateRGrid(IntVec<D> const & dimensions);

      /**
      * De-allocate fields in rgrid format.
      */
      void deallocateRGrid();

      /**
      * Allocate or re-allocate memory for fields in basis format.
      *
      * \param nBasis  number of basis functions 
      */
      void allocateBasis(int nBasis);

      /**
      * De-allocate fields in basis format.
      */
      void deallocateBasis();

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

      /**
      * Set field component values, in symmetrized Fourier format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
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
      * \param fields  array of new fields in r-grid format
      * \param isSymmetric is this field symmetric under the space group?
      */
      void setRGrid(DArray< RField<D> > const & fields, 
                    bool isSymmetric = false);

      /**
      * Read field component values from input stream, in symmetrized 
      * Fourier format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      * 
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      *
      * \param in  input stream from which to read fields
      * \param unitCell  associated crystallographic unit cell
      */
      void readBasis(std::istream& in, UnitCell<D>& unitCell);

      /**
      * Read field component values from file, in symmetrized 
      * Fourier format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      * 
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      *
      * \param filename  file from which to read fields
      * \param unitCell  associated crystallographic unit cell
      */
      void readBasis(std::string filename, UnitCell<D>& unitCell);

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
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      * 
      * \param in  input stream from which to read fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::istream& in, UnitCell<D>& unitCell,
                     bool isSymmetric = false);

      /**
      * Reads fields from a file in real-space (r-grid) format.
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
      * This object must already be allocated and associated with
      * a FieldIo object to run this function.
      * 
      * \param filename  file from which to read fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  is this field symmetric under the space group?
      */
      void readRGrid(std::string filename, UnitCell<D>& unitCell,
                     bool isSymmetric = false);

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
      DArray< RField<D> > const & rgrid() const;

      /**
      * Get the field for one monomer type in r-space grid format.
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      RField<D> const & rgrid(int monomerId) const;

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
      * This flat is set true if the fields were input in basis format 
      * by the function setBasis, or if they were set in grid format
      * by the function setRGrid but isSymmetric was set true.
      */
      bool isSymmetric() const;

   private:

      /*
      * Array of fields in symmetry-adapted basis format
      *
      * Element basis_[i] is an array that contains the components
      * of the field associated with monomer i, in a symmetry-adapted
      * Fourier expansion. 
      */
      DArray< DArray<double> > basis_;

      /*
      * Array of fields in real-space grid (r-grid) format
      *
      * Element basis_[i] is an RField<D> that contains values of the 
      * field associated with monomer i on the nodes of a regular mesh.
      */
      DArray< RField<D> > rgrid_;

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
      * Number of monomer types (number of fields)
      */
      int nMonomer_;

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

   };

   // Inline member functions

   // Get array of all fields in basis format (const)
   template <int D>
   inline
   DArray< DArray<double> > const & WFieldContainer<D>::basis() const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_;
   }

   // Get one field in basis format (const)
   template <int D>
   inline
   DArray<double> const & WFieldContainer<D>::basis(int id) const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_[id];
   }

   // Get all fields in r-grid format (const)
   template <int D>
   inline
   DArray< RField<D> > const &
   WFieldContainer<D>::rgrid() const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_;
   }

   // Get one field in r-grid format (const)
   template <int D>
   inline
   RField<D> const & WFieldContainer<D>::rgrid(int id) const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_[id];
   }

   // Has memory been allocated for fields in r-grid format?
   template <int D>
   inline bool WFieldContainer<D>::isAllocatedRGrid() const
   {  return isAllocatedRGrid_; }

   // Has memory been allocated for fields in basis format?
   template <int D>
   inline bool WFieldContainer<D>::isAllocatedBasis() const
   {  return isAllocatedBasis_; }

   // Has field data been initialized ?
   template <int D>
   inline bool WFieldContainer<D>::hasData() const
   {  return hasData_; }

   // Are the fields symmetric under space group operations?
   template <int D>
   inline bool WFieldContainer<D>::isSymmetric() const
   {  return isSymmetric_; }

   #ifndef PSPC_W_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class WFieldContainer<1>;
   extern template class WFieldContainer<2>;
   extern template class WFieldContainer<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
