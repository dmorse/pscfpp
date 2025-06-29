#ifndef RPC_BASIS_FIELD_STATE_TPP
#define RPC_BASIS_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"
#include "FieldState.tpp"
#include <util/global.h>

namespace Pscf {
namespace Rpc
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   BasisFieldState<D>::BasisFieldState()
    : FieldState<D, DArray<double> >()
   {}
  
   /*
   * Constructor.
   */
   template <int D>
   BasisFieldState<D>::BasisFieldState(System<D>& system)
    : FieldState<D, DArray<double> >(system)
   {}

   /*
   * Destructor.
   */
   template <int D>
   BasisFieldState<D>::~BasisFieldState()
   {}

   /*
   * Allocate all fields.
   */
   template <int D>
   void BasisFieldState<D>::allocate()
   {
      // Precondition
      UTIL_CHECK(hasSystem());

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
   template <int D>
   void BasisFieldState<D>::read(const std::string & filename)
   {
      allocate();
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      fieldIo.readFieldsBasis(filename, 
                                                    fields(), 
                                                    unitCell());
   }

   /**
   * Write fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::write(const std::string & filename)
   {
      FieldIo<D> const & fieldIo = system().domain().fieldIo();
      fieldIo.writeFieldsBasis(filename, fields(), unitCell());
   }

   /*
   * Gjt current state of associated System.
   */
   template <int D>
   void BasisFieldState<D>::getSystemState()
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
   template <int D>
   void BasisFieldState<D>::setSystemState(bool isFlexible)
   {
      system().w().setBasis(fields());
      if (isFlexible) {
         system().setUnitCell(unitCell());
      }
   }

} // namespace Rpc
} // namespace Pscf
#endif
