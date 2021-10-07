#ifndef PSPC_BASIS_FIELD_STATE_TPP
#define PSPC_BASIS_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"
#include "FieldState.tpp"

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
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

   /**
   * Read fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::read(std::string & filename)
   {
      allocate();
      fieldIo().readFieldsBasis(filename, fields());
   }

   /**
   * Write fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::write(std::string & filename)
   {
      fieldIo().writeFieldsBasis(filename, fields());
   }

   /*
   * Get current state of associated System.
   */
   template <int D>
   void BasisFieldState<D>::getSystemState()
   {
      // Get system unit cell
      unitCell() = system().unitCell();

      // Get system wFields
      allocate();
      int nMonomer = system().mixture().nMonomer();
      int nStar    = system().basis().nStar();
      int i, j;
      for (i = 0; i < nMonomer; ++i) {
         DArray<double>& stateField = field(i);
         const DArray<double>& systemField = system().wField(i);
         for (j = 0; j < nStar; ++j) {
            stateField[j] = systemField[j];
         }
      }

   }

   /*
   * Get current state of associated System.
   */
   template <int D>
   void BasisFieldState<D>::setSystemState()
   {
      // Update system unitCell
      system().unitCell() = unitCell();

      // Update system  wFields
      int nMonomer = system().mixture().nMonomer();
      int nStar = system().basis().nStar();
      int i, j;
      for (i = 0; i < nMonomer; ++i) {
         const DArray<double>& stateField = field(i);
         DArray<double>& systemField = system().wField(i);
         for (j = 0; j < nStar; ++j) {
            systemField[j] = stateField[j];
         }
      }

      // Update system wFieldsRgrid
      system().fieldIo().convertBasisToRGrid(system().wFields(),
                                             system().wFieldsRGrid());
   }

   /*
   * Allocate fields if necessary.
   */
   template <int D>
   void BasisFieldState<D>::allocate()
   {
      int nMonomer = system().mixture().nMonomer();
      int nStar    = system().basis().nStar();
      if (fields().isAllocated()) {
         UTIL_CHECK(fields().capacity() == nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK(field(i).isAllocated());
            UTIL_CHECK(field(i).capacity() == nStar);
         }
      } else {
         fields().allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            field(i).allocate(nStar);
         }
      }
   }

} // namespace Pspc
} // namespace Pscf
#endif
