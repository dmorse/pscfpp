#ifndef RPC_BASIS_FIELD_STATE_H
#define RPC_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldState.h"
#include <string>

namespace Pscf {
namespace Rpc
{

   using namespace Util;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   */
   template <int D>
   class BasisFieldState : public FieldState<D, DArray<double> >
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
      * \param filename name of input w-field file in symmetry-adapted format.
      */
      void read(const std::string & filename);
   
      /**
      * Write state to file.
      *
      * \param filename name of output file, in symmetry-adapated format.
      */
      void write(const std::string & filename);

      /**
      * Copy the current state of the associated system.
      *
      * Copy the fields and the unit cell.
      */
      void getSystemState();

      /**
      * Set the state of the associated system to this state.
      *
      * \param isFlexible update system unit cell iff isFlexible == true.
      */
      void setSystemState(bool isFlexible);

      // Inherited member functions
      using FieldState<D, DArray<double> >::fields;
      using FieldState<D, DArray<double> >::field;
      using FieldState<D, DArray<double> >::unitCell;
      using FieldState<D, DArray<double> >::system;
      using FieldState<D, DArray<double> >::hasSystem;
      using FieldState<D, DArray<double> >::setSystem;

   };

   // Explicit instantiation declarations
   extern template class BasisFieldState<1>;
   extern template class BasisFieldState<2>;
   extern template class BasisFieldState<3>;

} // namespace Rpc
} // namespace Pscf
#endif
