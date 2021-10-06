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
   *    - a UnitCell
   *    - Monomer chemical fields in both basis and grid formats
   *
   * \ingroup Pscf_Pspc_Module
   */
   template <int D, class FT>
   class FieldState 
   {

   public:

      /// Field type
      typedef FT Field;

      /// \name Construction and Destruction
      //@{

      /**
      * Default constructor.
      */
      FieldState();

      /**
      * Constructor, creates association with a System.
      */
      FieldState(const System& system);

      /**
      * Destructor.
      */
      ~FieldState();

      /**
      * Get array of all chemical potential fields.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<Field>& fields();

      /**
      * Get field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      DArray<Field>& field(int monomerId);

      /**
      * Get UnitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D>& unitCell();

   protected:

      /**
      * Get associated System.
      */
      System<D>& system();

      /**
      * Get associated FieldIo
      */
      FieldIo<D>& fieldIo();

   private:

      /**
      * Array of fields for all monomer types.
      */
      DArray<Field> fields_;

      /**
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

      /**
      * Pointer to associated System.
      */
      System<D>* systemPtr_;

   };

   // Inline member functions

   template <int D, class FT>
   inline
   DArray<Field>& FieldState<D,FT>::fields()
   {  return fields_; }

   template <int D, class FT>
   inline
   DArray<Field>& FieldState<D,FT>::field(int id)
   {  return fields_[id]; }

   // Get the internal UnitCell<D> object.
   template <int D, class FT>
   inline UnitCell<D>& FieldState<D,FT>::unitCell()
   { return unitCell_; }

   // Get the internal FieldIo<D> object.
   template <int D, class FT>
   inline FieldIo<D>& FieldState<D,FT>::fieldIo()
   { return fieldIo_; }

   // Get the associated System<D> object.
   template <int D, class FT>
   inline System<D>& FieldState<D,FT>::system()
   { return *systemPtr_; }

   /**
   * FieldState for fields in symmetry-adapted basis format.
   */
   template <int D>
   class BasisFieldState : public FieldState<D, DArray<double> >
   {

      /**
      * Read state from file.
      *
      * \param filename name of input w-field basis file
      */
      void read(std::string & filename);
   
      /**
      * Write state to file.
      *
      * \param filename name of input w-field r-grid file
      */
      void write(std::string & filename);
   
   }

   #if 0
   #ifndef PSPC_FIELD_STATE_TPP
   // Suppress implicit instantiation
   extern template class BasisFieldState<1>;
   extern template class BasisFieldState<2>;
   extern template class BasisFieldState<3>;
   #endif
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
