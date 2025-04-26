#ifndef RPG_MASK_TPP
#define RPG_MASK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
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
   * Calculate the average value of the rgrid_ member.
   */
   template <int D>
   double Mask<D>::rGridAverage() const
   {
      RField<D> const & rg = MaskTmpl< D, FieldIo<D>, RField<D> >::rgrid();
      return (Reduce::sum(rg) / ((double)rg.capacity()));
   }

} // namespace Rpg
} // namespace Pscf
#endif
