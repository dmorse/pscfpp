#ifndef RP_BASIS_FIELD_STATE_TPP
#define RP_BASIS_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>
#include <prdc/crystal/Basis.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Default constructor.
   */
   template <int D, class T>
   BasisFieldState<D,T>::BasisFieldState()
    : FieldStateT()
   {}

   /*
   * Constructor.
   */
   template <int D, class T>
   BasisFieldState<D,T>::BasisFieldState(System<D,T>& system)
    : FieldStateT(system)
   {}

   /*
   * Allocate all fields.
   */
   template <int D, class T>
   void BasisFieldState<D,T>::allocate()
   {
      // Precondition
      UTIL_CHECK(FieldStateT::hasSystem());

      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      if (fields().isAllocated()) {
         UTIL_CHECK(fields().capacity() == nMonomer);
      } else {
         fields().allocate(nMonomer);
      }

      int nBasis = system().domain().basis().nBasis();
      UTIL_CHECK(nBasis > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (field(i).isAllocated()) {
            UTIL_CHECK(field(i).capacity() == nBasis);
         } else {
            field(i).allocate(nBasis);
         }
      }

   }

   /*
   * Read fields in symmetry-adapted basis format.
   */
   template <int D, class T>
   void BasisFieldState<D,T>::read(const std::string & filename)
   {
      allocate();
      FieldIoT const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsBasis(filename, fields(), unitCell());
   }

   /*
   * Write fields in symmetry-adapted basis format.
   */
   template <int D, class T>
   void BasisFieldState<D,T>::write(const std::string & filename)
   {
      FieldIoT const & fieldIo = system().domain().fieldIo();
      fieldIo.writeFieldsBasis(filename, fields(), unitCell());
   }

   /*
   * Get current state of associated System.
   */
   template <int D, class T>
   void BasisFieldState<D,T>::getSystemState()
   {
      // Get system wFields
      allocate();
      int nMonomer = system().mixture().nMonomer();
      int nBasis = system().domain().basis().nBasis();
      int i, j;
      for (i = 0; i < nMonomer; ++i) {
         DArray<double>& stateField = field(i);
         const DArray<double>& systemField = system().w().basis(i);
         for (j = 0; j < nBasis; ++j) {
            stateField[j] = systemField[j];
         }
      }

      // Get system unit cell
      unitCell() = system().domain().unitCell();
   }

   /*
   * Set System state to current state of the BasisFieldState object.
   */
   template <int D, class T>
   void BasisFieldState<D,T>::setSystemState(bool newCellParams)
   {
      system().w().setBasis(fields());
      if (newCellParams) {
         system().setUnitCell(unitCell());
      }
   }

} // namespace Rp
} // namespace Pscf
#endif
