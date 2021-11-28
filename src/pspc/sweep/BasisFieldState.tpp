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

      int nStar = system().basis().nStar();
      UTIL_CHECK(nStar > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (field(i).isAllocated()) {
            UTIL_CHECK(field(i).capacity() == nStar);
         } else {
            field(i).allocate(nStar);
         }
      }

   }
 
   /**
   * Read fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::read(const std::string & filename)
   {
      allocate();
      system().fieldIo().readFieldsBasis(filename, fields(), unitCell());
   }

   /**
   * Write fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::write(const std::string & filename)
   {
      system().fieldIo().writeFieldsBasis(filename, fields(), unitCell());
   }

   /*
   * Get current state of associated System.
   */
   template <int D>
   void BasisFieldState<D>::getSystemState()
   {
      std::cout << "Checkpoint 1.2! \n";
      // Get system unit cell
      unitCell() = system().unitCell();
      std::cout << "Checkpoint 1.3! \n";
      // Get system wFields
      allocate();
      int nMonomer = system().mixture().nMonomer();
      int nStar    = system().basis().nStar();
      int i, j;
      std::cout << "Checkpoint 1.4! \n";
      for (i = 0; i < nMonomer; ++i) {
         DArray<double>& stateField = field(i);
         const DArray<double>& systemField = system().wField(i);
         for (j = 0; j < nStar; ++j) {
            stateField[j] = systemField[j];
         }
      }
      std::cout << "Checkpoint 1.5! \n";

   }

   /*
   * Set System state to current state of the BasisFieldState object.
   */
   template <int D>
   void BasisFieldState<D>::setSystemState(bool isFlexible)
   {
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

      if (isFlexible) {
         // Update system unitCell
         system().unitCell() = unitCell();
      }

   }

} // namespace Pspc
} // namespace Pscf
#endif
