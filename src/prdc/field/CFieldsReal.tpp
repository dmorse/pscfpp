#ifndef PRDC_C_FIELDS_REAL_TPP
#define PRDC_C_FIELDS_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFieldsReal.h"
#include <prdc/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class RFT, class FIT>
   CFieldsReal<D,RFT,FIT>::CFieldsReal()
    : basis_(),
      rgrid_(),
      nMonomer_(0),
      writeUnitCellPtr_(nullptr),
      fieldIoPtr_(nullptr),
      isAllocatedRGrid_(false),
      isAllocatedBasis_(false)
   {}

   /*
   * Destructor.
   */
   template <int D, class RFT, class FIT>
   CFieldsReal<D,RFT,FIT>::~CFieldsReal()
   {}

   /*
   * Create an association with a FIT object.
   */
   template <int D, class RFT, class FIT>
   void
   CFieldsReal<D,RFT,FIT>::setFieldIo(FIT const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Set the stored value of nMonomer (this may only be called once).
   */
   template <int D, class RFT, class FIT>
   void CFieldsReal<D,RFT,FIT>::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
   }

   /*
   * Set the unit cell for parameters written to a field header.
   */
   template <int D, class RFT, class FIT>
   void CFieldsReal<D,RFT,FIT>::setWriteUnitCell(UnitCell<D> const & cell)
   {
      UTIL_CHECK(!writeUnitCellPtr_);
      writeUnitCellPtr_ = &cell;
   }

   /*
   * Allocate memory for fields in r-grid format.
   */
   template <int D, class RFT, class FIT>
   void
   CFieldsReal<D,RFT,FIT>::allocateRGrid(IntVec<D> const& meshDimensions)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(!isAllocatedRGrid_);

      // Allocate arrays
      rgrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].allocate(meshDimensions);
      }
      isAllocatedRGrid_ = true;
   }

   /*
   * Allocate memory for fields in basis format.
   */
   template <int D, class RFT, class FIT>
   void CFieldsReal<D,RFT,FIT>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(!isAllocatedBasis_);

      // Allocate
      basis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;
   }

   /*
   * Allocate memory for all fields.
   */
   template <int D, class RFT, class FIT>
   void 
   CFieldsReal<D,RFT,FIT>::allocate(int nMonomer, int nBasis,
                                       IntVec<D> const & meshDimensions)
   {
      setNMonomer(nMonomer);
      allocateRGrid(meshDimensions);
      allocateBasis(nBasis);
   }

} // namespace Prdc
} // namespace Pscf
#endif
