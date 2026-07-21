/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"                  // class header
#include <rpg/field/FieldIo.h>
#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/DeviceArray.h>
#include <rp/field/WFieldsBase.tpp>   // base class implementation

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   // Public member function

   /*
   * Set new w-field values, using unfolded array of r-grid fields.
   */
   template <int D>
   void WFields<D, CudaTp<D> >::setRGrid(DeviceArray<cudaReal>& fields)
   {
      // Create DArray tmp with RField<D, CudaTp<D> > elements
      DArray< RField<D, CudaTp<D> > > tmp;
      const int nMonomer = RpWFields::nMonomer();
      tmp.allocate(nMonomer);

      // Associate each RField<D, CudaTp<D> > with a slice of the unfolded array
      IntVec<D> const & meshDimensions = RpWFields::meshDimensions();
      const int meshSize = RpWFields::meshSize();
      for (int i = 0; i < nMonomer; i++) {
         tmp[i].associate(fields, i*meshSize, meshDimensions);
      }

      // Use tmp array to set w-fields for all monomer types
      bool isSymmetric = false;
      RpWFields::setRGrid(tmp, isSymmetric);
   }

} // namespace Rp
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class WFieldsBase<1, CudaTp<1> >;
      template class WFieldsBase<2, CudaTp<2> >;
      template class WFieldsBase<3, CudaTp<3> >;
      template class WFields<1, CudaTp<1> >;
      template class WFields<2, CudaTp<2> >;
      template class WFields<3, CudaTp<3> >;
   }
}
