#ifndef PRDC_FIELD_STATE_H
#define PRDC_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>        // member
#include <util/containers/DArray.h>       // member template

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Record of a state of a periodic system (fields + unit cell).
   * 
   * Template parameters:
   *   
   *    - D : dimensions of space
   *    - FT : a field or array type
   *    - ST : a system type
   *
   * A FieldState<D, FT, ST> has:
   *
   *    - An array of field objects of class FT
   *    - a UnitCell<D> object
   *    - a pointer to a system of class ST
   * 
   * Specializations of FieldState can be used to store either chemical 
   * potential or concentration fields, along with an associated unit
   * cell. Different choices for class FT can be used to store fields 
   * in symmetry-adapted basis function, r-grid or k-grid format.
   *
   * FieldState is a standard class template, in which all member
   * function definitions are located in this header file.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class FT, class ST>
   class FieldState 
   {

   public:

      /// \name Construction and Destruction
      ///@{

      /**
      * Default constructor.
      */
      FieldState();

      /**
      * Constructor, creates association with a System.
      *
      * Equivalent to default construction followed by setSystem(system).
      *
      * \param system associated parent ST object.
      */
      FieldState(ST& system);

      /**
      * Destructor.
      */
      ~FieldState();

      /**
      * Set association with System, after default construction.
      *
      * \param system associated parent ST object.
      */
      void setSystem(ST& system);

      ///@}
      /// \name Accessors
      ///@{

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

      ///@}

   protected:

      /**
      * Has a system been set?
      */
      bool hasSystem();

      /**
      * Get associated System by reference.
      */
      ST& system();

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
      ST* systemPtr_;

   };

   // Public inline member functions

   // Get an array of all fields (const reference)
   template <int D, class FT, class ST> inline
   const DArray<FT>& FieldState<D,FT,ST>::fields() const
   {  return fields_; }

   // Get an array of all fields (non-const reference)
   template <int D, class FT, class ST> inline
   DArray<FT>& FieldState<D,FT,ST>::fields()
   {  return fields_; }

   // Get field for monomer type id (const reference)
   template <int D, class FT, class ST> inline
   const FT& FieldState<D,FT,ST>::field(int id) const
   {  return fields_[id]; }

   // Get field for monomer type id (non-const reference)
   template <int D, class FT, class ST> inline 
   FT& FieldState<D,FT,ST>::field(int id)
   {  return fields_[id]; }

   // Get the internal Unitcell (const reference)
   template <int D, class FT, class ST> inline 
   const UnitCell<D>& FieldState<D,FT,ST>::unitCell() const
   { return unitCell_; }

   // Get the internal Unitcell (non-const reference)
   template <int D, class FT, class ST> inline 
   UnitCell<D>& FieldState<D,FT,ST>::unitCell()
   { return unitCell_; }

   // Protected inline member functions

   // Has the system been set?
   template <int D, class FT, class ST> inline 
   bool FieldState<D,FT,ST>::hasSystem()
   {  return (systemPtr_ != 0); }
   
   // Get the associated ST object.
   template <int D, class FT, class ST> inline 
   ST& FieldState<D,FT,ST>::system()
   {
      assert(systemPtr_);  
      return *systemPtr_; 
   }

   // Noninline public member functions

   /*
   * Constructor.
   */
   template <int D, class FT, class ST>
   FieldState<D,FT,ST>::FieldState()
    : fields_(),
      unitCell_(),
      systemPtr_(nullptr)
   {}
  
   /*
   * Constructor.
   */
   template <int D, class FT, class ST>
   FieldState<D,FT,ST>::FieldState(ST& system)
    : fields_(),
      unitCell_(),
      systemPtr_(nullptr)
   {  setSystem(system); }

   /*
   * Destructor.
   */
   template <int D, class FT, class ST>
   FieldState<D,FT,ST>::~FieldState()
   {}

   /*
   * Set association with system, after default construction.
   */
   template <int D, class FT, class ST>
   void FieldState<D,FT,ST>::setSystem(ST& system)
   {
      if (hasSystem()) {
         UTIL_CHECK(systemPtr_ == &system);
      } else {
         systemPtr_ = &system;
      }
   }

} // namespace Prdc
} // namespace Pscf
#endif
