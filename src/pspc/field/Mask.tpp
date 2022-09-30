#ifndef PSPC_MASK_TPP
#define PSPC_MASK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <pspc/field/FieldIo.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Mask<D>::Mask()
    : basis_(),
      rgrid_(),
      fieldIoPtr_(0),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      isAllocated_(false),
      hasData_(false),
      isSymmetric_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   Mask<D>::~Mask()
   {}

   /*
   * Create an association with a FieldIo object.
   */
   template <int D>
   void Mask<D>::setFieldIo(FieldIo<D> const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Allocate memory for field.
   */
   template <int D>
   void Mask<D>::allocate(int nBasis, IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(!isAllocated_);

      // Set mesh and basis dimensions
      nBasis_ = nBasis;
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }
  
      // Allocate field arrays 
      basis_.allocate(nBasis);
      rgrid_.allocate(meshDimensions);
      isAllocated_ = true;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void Mask<D>::setBasis(DArray<double> const & field)
   {
      UTIL_CHECK(field.capacity() == nBasis_);
      for (int j = 0; j < nBasis_; ++j) {
         basis_[j] = field[j];
      }
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);
      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Set new field values, using r-grid field as inputs.
   */
   template <int D>
   void Mask<D>::setRGrid(RField<D> const & field,
                          bool isSymmetric)
   {
      for (int j = 0; j < meshSize_; ++j) {
         rgrid_[j] = field[j];
      }
      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }
      hasData_ = true;
      isSymmetric_ =  isSymmetric;
   }

} // namespace Pspc
} // namespace Pscf
#endif
