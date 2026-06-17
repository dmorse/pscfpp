/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"                  // class header
#include <rpg/field/FieldIo.h>
#include <prdc/cuda/RField.h>
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
      // Create DArray tmp with RField<D> elements
      DArray< RField<D> > tmp;
      const int nMonomer = RpWFields::nMonomer();
      tmp.allocate(nMonomer);

      // Associate each RField<D> with a slice of the unfolded array
      IntVec<D> const & meshDimensions = RpWFields::meshDimensions();
      const int meshSize = RpWFields::meshSize();
      for (int i = 0; i < nMonomer; i++) {
         tmp[i].associate(fields, i*meshSize, meshDimensions);
      }

      // Use tmp array to set w-fields for all monomer types
      bool isSymmetric = false;
      RpWFields::setRGrid(tmp, isSymmetric);
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      template class WFields<1, Rpg::Types<1> >;
      template class WFields<2, Rpg::Types<2> >;
      template class WFields<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class WFields<1>;
      template class WFields<2>;
      template class WFields<3>;
   }
}
