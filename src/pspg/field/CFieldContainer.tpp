#ifndef PSPG_C_FIELD_CONTAINER_TPP
#define PSPG_C_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFieldContainer.h"

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   CFieldContainer<D>::CFieldContainer()
    : basis_(),
      rgrid_(),
      nMonomer_(0),
      isAllocatedRGrid_(false),
      isAllocatedBasis_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   CFieldContainer<D>::~CFieldContainer()
   {}

   /*
   * Set the stored value of nMonomer (this may only be called once).
   */
   template <int D>
   void CFieldContainer<D>::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
   }

   /*
   * Allocate memory for fields in r-grid format.
   */
   template <int D>
   void
   CFieldContainer<D>::allocateRGrid(IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(nMonomer_ > 0);

      // If already allocated, deallocate.
      if (isAllocatedRGrid_) {
         deallocateRGrid();
      }

      // Allocate arrays
      rgrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].allocate(meshDimensions);
      }
      isAllocatedRGrid_ = true;
   }

   /*
   * De-allocate memory for fields in r-grid format
   */
   template <int D>
   void CFieldContainer<D>::deallocateRGrid()
   {
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(nMonomer_ > 0);
      for (int i = 0; i < nMonomer_ ; ++i) {
         rgrid_[i].deallocate();
      }
      rgrid_.deallocate();
      isAllocatedRGrid_ = false;
   }

   /*
   * Allocate memory for fields in basis format.
   */
   template <int D>
   void CFieldContainer<D>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(nMonomer_ > 0);

      // If already allocated, deallocate.
      if (isAllocatedBasis_) {
         deallocateBasis();
      }

      // Allocate
      basis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;
   }

   /*
   * De-allocate memory for fields in basis format.
   */
   template <int D>
   void CFieldContainer<D>::deallocateBasis()
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(isAllocatedBasis_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].deallocate();
      }
      basis_.deallocate();
      isAllocatedBasis_ = false;
   }

   /*
   * Allocate memory for all fields.
   */
   template <int D>
   void CFieldContainer<D>::allocate(int nMonomer, int nBasis,
                                    IntVec<D> const & meshDimensions)
   {
      setNMonomer(nMonomer);
      allocateRGrid(meshDimensions);
      allocateBasis(nBasis);
   }

} // namespace Pspg
} // namespace Pscf
#endif
