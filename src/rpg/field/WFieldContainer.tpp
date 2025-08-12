#ifndef RPG_W_FIELD_CONTAINER_TPP
#define RPG_W_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldContainer.h"
#include <prdc/field/WFieldsReal.tpp>
#include <prdc/cuda/VecOp.h>
#include <pscf/cuda/DeviceArray.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   // Public member function

   /*
   * Set new w-field values, using unfolded array of r-grid fields.
   */
   template <int D>
   void WFieldContainer<D>::setRGrid(DeviceArray<cudaReal>& fields)
   {
      DArray< RField<D> > tmp;
      tmp.allocate(nMonomer());
      for (int i = 0; i < nMonomer(); i++) {
         tmp[i].associate(fields, i * meshSize(), meshDimensions());
      }
      bool isSymmetric = false;
      Base::setRGrid(tmp, isSymmetric);
   }

   // Private virtual function

   template <int D>
   void 
   WFieldContainer<D>::assignRField(RField<D>& lhs, RField<D> const & rhs) 
   const
   {
      int n = rhs.capacity();
      UTIL_CHECK(lhs.capacity() == n);
      UTIL_CHECK(meshSize() == n);
      VecOp::eqV(lhs, rhs);
   }

} // namespace Rpg
} // namespace Pscf
#endif
