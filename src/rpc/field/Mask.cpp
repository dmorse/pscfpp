/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <pscf/cpu/Reduce.h>

#include <rp/field/Mask.tpp>    // base class template implementation

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /*
   * Return volume fraction of the unit cell occupied by the
   * polymers/solvents.
   */
   template <int D>
   double Mask<D>::rGridAverage() const
   {
      RField<D> const & rg = Rp::Mask< D, RField<D>, FieldIo<D> >::rgrid();

      // Sum up elements of r-grid mask field.
      // Use Kahan summation to reduce accumulation of error
      double sum(0.0), err(0.0), tempVal, tempSum;
      int n = rg.capacity();
      for (int i = 0; i < n; ++i) {
         tempVal = rg[i] - err;
         tempSum = sum + tempVal;
         err = tempSum - sum - tempVal;
         sum = tempSum;
      }

      return (sum / ((double)rg.capacity()));
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cpu;
      template class Rp::Mask< 1, RField<1>, Rpc::FieldIo<1> >;
      template class Rp::Mask< 2, RField<2>, Rpc::FieldIo<2> >;
      template class Rp::Mask< 3, RField<3>, Rpc::FieldIo<3> >;
   }
   namespace Rpc {
      template class Mask<1>;
      template class Mask<2>;
      template class Mask<3>;
   }
}
