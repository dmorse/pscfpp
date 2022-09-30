#ifndef PSPC_C_FIELD_CONTAINER_TPP
#define PSPC_C_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFieldContainer.h"

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   CFieldContainer<D>::CFieldContainer()
    : basis_(),
      rgrid_(),
      isAllocated_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   CFieldContainer<D>::~CFieldContainer()
   {}

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void CFieldContainer<D>::allocate(int nMonomer, int nBasis, 
                                    IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(!isAllocated_);
      basis_.allocate(nMonomer);
      rgrid_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         basis_[i].allocate(nBasis);
         rgrid_[i].allocate(meshDimensions);
      }
      isAllocated_ = true;
   }

} // namespace Pspc
} // namespace Pscf
#endif
