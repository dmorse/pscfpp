#ifndef RPG_MASK_TPP
#define RPG_MASK_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.h"
#include <prdc/cuda/Reduce.h>

namespace Pscf {
namespace Rpg
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Calculate the average value of the rgrid_ member.
   */
   template <int D>
   double Mask<D>::rGridAverage() const
   {
      RField<D> const & rg = MaskReal< D, RField<D>, FieldIo<D> >::rgrid();
      return (Reduce::sum(rg) / ((double)rg.capacity()));
   }

} // namespace Rpg
} // namespace Pscf
#endif
