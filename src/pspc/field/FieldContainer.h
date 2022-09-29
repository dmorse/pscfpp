#ifndef PSPC_FIELD_CONTAINER_H
#define PSPC_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>     // base class

#include <pspc/field/RField.h>             // member template parameter
#include <pscf/math/IntVec.h>              // function parameter
#include <util/containers/DArray.h>        // member template

namespace Pscf {
namespace Pspc {

   template <int D> class FieldIo;

   using namespace Util;

   /**
   * A list of fields stored in both basis and r-grid format.
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class FieldContainer : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      FieldContainer();

      /**
      * Destructor.
      */
      ~FieldContainer();

      /**
      * Create association with FieldIo (store pointer).
      */
      void setFieldIo(FieldIo<D> const & fieldIo);

      /**
      * Allocate memory for fields.
      */
      void allocate(int nMonomer, int nBasis, IntVec<D> const & dimensions);

      /**
      * Set field component values, in symmetrized Fourier format.
      *
      * This function also computes and stores the corresponding
      * r-grid representation. On return, hasData and isSymmetric
      * are both true.
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setBasis(DArray< DArray<double> > const & fields);

      /**
      * Set fields values in real-space (r-grid) format.
      *
      * If the isSymmetric parameter is true, this function assumes that 
      * the fields are known to be symmetric and so computes and stores
      * the corresponding basis components.
      * 
      * \param fields  array of new w (chemical potential) fields
      * \param isSymmetric if true, attempt to compute basis format
      */
      void setRGrid(DArray< RField<D> > const & fields, 
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
      * \param monomerId integer monomer type index
      */
      DArray<double> const & basis(int monomerId) const;

      /**
      * Get array of all fields in r-space grid format.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< RField<D> > const & rgrid() const;

      /**
      * Get the field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      RField<D> const & rgrid(int monomerId) const;

      /**
      * Have the fields been set?
      *
      * This flag is set true in setBasis and setRGrid.
      */
      bool hasData() const;

      /**
      * Arefields symmetric under all elements of the space group?
      *
      * This is set true if the fields were input in basis format 
      * by the function setBasis, or if they were set in grid format
      * by the function setRGrid but isSymmetric was set true.
      */
      bool isSymmetric() const;

   private:

      /*
      * Array of fields in symmetry-adapted basis format
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< DArray<double> > basis_;

      /*
      * Array of fields in real-space grid (r-grid) format
      *
      * Indexed by monomer typeId, size = nMonomer.
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
      * Has memory been allocated for fields?
      */
      bool isAllocated_;

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

   // Get array of all monomer chemical potential fields.
   template <int D>
   inline
   DArray< DArray<double> > const & FieldContainer<D>::basis() const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_;
   }

   // Get one monomer chemical potential field.
   template <int D>
   inline
   DArray<double> const & FieldContainer<D>::basis(int id) const
   {
      UTIL_ASSERT(hasData_);
      UTIL_ASSERT(isSymmetric_);
      return basis_[id];
   }

   // Get an array of monomer chemical potential fields on r-space grids.
   template <int D>
   inline
   DArray< RField<D> > const &
   FieldContainer<D>::rgrid() const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_;
   }

   // Get a single monomer chemical potential field on an r-space grid.
   template <int D>
   inline
   RField<D> const & FieldContainer<D>::rgrid(int id) const
   {
      UTIL_ASSERT(hasData_);
      return rgrid_[id];
   }

   // Have the field data been set?
   template <int D>
   inline bool FieldContainer<D>::hasData() const
   {  return hasData_; }

   // Are the fields symmetric under space group operations?
   template <int D>
   inline bool FieldContainer<D>::isSymmetric() const
   {  return isSymmetric_; }

   #ifndef PSPC_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class FieldContainer<1>;
   extern template class FieldContainer<2>;
   extern template class FieldContainer<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
