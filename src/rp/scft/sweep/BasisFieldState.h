#ifndef RP_BASIS_FIELD_STATE_H
#define RP_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/fieldIo/FieldState.h>      // base class
#include <pscf/backends/TmplDeclare.h>    // preprocessor macros

#include <string>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class FieldIo;
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   *
   * Template parameters:
   *
   *   - D : dimension of space (1, 2, or 3)
   *   - T : backend identifier class (e.g. (CPT or CUT))
   *
   * \ingroup Rp_Scft_Sweep_Module
   */
   template <int D, class T>
   class BasisFieldState
    : public FieldState<D, DArray<double>, System<D,T> >
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
      BasisFieldState(System<D,T>& system);

      ~BasisFieldState() = default;

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

      using FieldStateT = FieldState<D, DArray<double>, System<D,T> >;

      // Inherited member functions
      using FieldStateT::fields;
      using FieldStateT::field;
      using FieldStateT::unitCell;
      using FieldStateT::system;

   private:

      using FieldIoT = FieldIo<D,T>;

   };

   // Explicit instantation declarations
   PSCF_TMPL_DECLARE(BasisFieldState)

} // namespace Rp
} // namespace Pscf

#endif
