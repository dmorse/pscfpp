#ifndef RPG_BASIS_FIELD_STATE_H
#define RPG_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/FieldState.h>
#include <string>

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Prdc;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   */
   template <int D>
   class BasisFieldState
    : public FieldState<D, DArray<double>, System<D> >
   {
   public:

      /**
      * Default constructor.
      */
      BasisFieldState();

      /**
      * Constructor, create association with a parent system.
      *
      * \param system associated parent system
      */
      BasisFieldState(System<D>& system);

      /**
      * Destructor.
      */
      ~BasisFieldState();

      /**
      * Allocate all fields.
      *
      * Precondition: hasSystem() == true
      */
      void allocate();

      /**
      * Read state from file.
      *
      * \param filename  name of input w-field file, in basis format
      */
      void read(const std::string & filename);

      /**
      * Write state to file.
      *
      * \param filename  name of output w-field file, in basis format
      */
      void write(const std::string & filename);

      /**
      * Store the current state of the associated system.
      *
      * Copy the fields and the unit cell.
      */
      void getSystemState();

      /**
      * Set the state of the associated system to this state.
      *
      * \param newCellParams  update unit cell parameters iff true
      */
      void setSystemState(bool newCellParams);

      // Inherited member functions
      using FieldStateT = FieldState<D, DArray<double>, System<D> >;
      using FieldStateT::fields;
      using FieldStateT::field;
      using FieldStateT::unitCell;
      using FieldStateT::system;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Prdc {
      extern template
      class FieldState< 1, DArray<double>, Rpg::System<1> >;
      extern template
      class FieldState< 2, DArray<double>, Rpg::System<2> >;
      extern template
      class FieldState< 3, DArray<double>, Rpg::System<3> >;
   }
   namespace Rpg {
      extern template class BasisFieldState<1>;
      extern template class BasisFieldState<2>;
      extern template class BasisFieldState<3>;
   }
}
#endif
