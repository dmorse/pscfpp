#ifndef PSPC_FIELD_STATE_H
#define PSPC_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/UnitCell.h>         // member
#include <pspc/field/FieldIo.h>            // member
#include <util/containers/DArray.h>        // member template

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   template <int D> class System;

   /**
   * Record of a state of a System (fields + unit cell).
   *
   *    - a UnitCell<D> object
   *    - An array of field objects of class FT
   *
   * The template parameter D is the dimension of space, while
   * parameter FT is a field type.
   *
   * A FieldState can be used to store either chemical potential or
   * concentration fields, along with an associated UnitCell<D>. 
   * Different choices for class FT can be used to store fields in
   * symmetry-adapted basis function, r-grid or k-grid format.
   *
   * \ingroup Pspc_Sweep_Module
   */
   template <int D, class FT>
   class FieldState 
   {

   public:

      /// \name Construction and Destruction
      //@{

      /**
      * Default constructor.
      */
      FieldState();

      /**
      * Constructor, creates association with a System.
      *
      * Equivalent to default construction followed by setSystem(system).
      *
      * \param system associated parent System<D> object.
      */
      FieldState(System<D>& system);

      /**
      * Destructor.
      */
      ~FieldState();

      /**
      * Set association with System, after default construction.
      *
      * \param system associated parent System<D> object.
      */
      void setSystem(System<D>& system);

      //@}
      /// \name Accessors
      //@{

      /**
      * Get array of all fields by const reference.
      *
      * The array capacity is equal to the number of monomer types.
      */
      const DArray<FT>& fields() const;

      /**
      * Get array of all chemical potential fields (non-const reference).
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<FT>& fields();

      /**
      * Get a field for a single monomer type by const reference.
      *
      * \param monomerId integer monomer type index
      */
      const FT& field(int monomerId) const;

      /**
      * Get field for a specific monomer type (non-const reference).
      *
      * \param monomerId integer monomer type index
      */
      FT& field(int monomerId);

      /**
      * Get UnitCell (i.e., lattice type and parameters) by const reference.
      */
      const UnitCell<D>& unitCell() const;

      /**
      * Get the UnitCell by non-const reference.
      */
      UnitCell<D>& unitCell();

      //@}

   protected:

      /**
      * Has a system been set?
      */
      bool hasSystem();

      /**
      * Get associated System by reference.
      */
      System<D>& system();

   private:

      /**
      * Array of fields for all monomer types.
      */
      DArray<FT> fields_;

      /**
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * Pointer to associated system.
      */
      System<D>* systemPtr_;

   };

   // Public inline member functions

   // Get an array of all fields (const reference)
   template <int D, class FT>
   inline
   const DArray<FT>& FieldState<D,FT>::fields() const
   {  return fields_; }

   // Get an array of all fields (non-const reference)
   template <int D, class FT>
   inline
   DArray<FT>& FieldState<D,FT>::fields()
   {  return fields_; }

   // Get field for monomer type id (const reference)
   template <int D, class FT>
   inline
   const FT& FieldState<D,FT>::field(int id) const
   {  return fields_[id]; }

   // Get field for monomer type id (non-const reference)
   template <int D, class FT>
   inline FT& FieldState<D,FT>::field(int id)
   {  return fields_[id]; }

   // Get the internal Unitcell (const reference)
   template <int D, class FT>
   inline 
   const UnitCell<D>& FieldState<D,FT>::unitCell() const
   { return unitCell_; }

   // Get the internal Unitcell (non-const reference)
   template <int D, class FT>
   inline UnitCell<D>& FieldState<D,FT>::unitCell()
   { return unitCell_; }

   // Protected inline member functions

   // Has the system been set?
   template <int D, class FT>
   inline bool FieldState<D,FT>::hasSystem()
   {  return (systemPtr_ != 0); }
   
   // Get the associated System<D> object.
   template <int D, class FT>
   inline System<D>& FieldState<D,FT>::system()
   {
      assert(systemPtr_ != 0);  
      return *systemPtr_; 
   }

   #ifndef PSPC_FIELD_STATE_TPP
   // Suppress implicit instantiation
   extern template class FieldState< 1, DArray<double> >;
   extern template class FieldState< 2, DArray<double> >;
   extern template class FieldState< 3, DArray<double> >;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
