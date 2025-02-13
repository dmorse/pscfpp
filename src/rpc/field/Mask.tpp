#ifndef RPC_MASK_TPP
#define RPC_MASK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <prdc/cpu/RFieldDft.h>

namespace Pscf {
namespace Rpc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Mask<D>::Mask()
    : MaskTmpl< D, FieldIo<D>, RField<D> >()
   {}

   /*
   * Destructor.
   */
   template <int D>
   Mask<D>::~Mask()
   {}

   /*
   * Return volume fraction of the unit cell occupied by the 
   * polymers/solvents.
   */
   template <int D>
   double Mask<D>::rGridAverage() const
   {
      RField<D> const & rg = MaskTmpl< D, FieldIo<D>, RField<D> >::rgrid();

      // Sum up elements of rg.
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

} // namespace Rpc
} // namespace Pscf
#endif
