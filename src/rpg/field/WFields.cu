/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"                  // class header
#include <rpg/field/FieldIo.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/DeviceArray.h>
#include <rp/field/WFields.tpp>       // base class implementation

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
   void WFields<D>::setRGrid(DeviceArray<cudaReal>& fields)
   {
      int nMonomer = Base::nMonomer();
      int meshSize = Base::meshSize();
      DArray< RField<D> > tmp;
      tmp.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         tmp[i].associate(fields, i * meshSize, Base::meshDimensions());
      }

      bool isSymmetric = false;
      Base::setRGrid(tmp, isSymmetric);
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      template class WFields<1, Cuda::RField<1>, Rpg::FieldIo<1> >;
      template class WFields<2, Cuda::RField<2>, Rpg::FieldIo<2> >;
      template class WFields<3, Cuda::RField<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      template class WFields<1>;
      template class WFields<2>;
      template class WFields<3>;
   } 
} 
