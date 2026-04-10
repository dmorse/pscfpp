#ifndef RPG_MASK_H
#define RPG_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Mask.h>       // base class template
#include "FieldIo.h"             // base class template argument
#include <prdc/cuda/RField.h>    // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * A field to which the total monomer concentration is constrained.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Mask, and
   * inherit their public interface and all of their source code from
   * this base class.  
   *
   * \see Rp::Mask
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Mask
     : public Rp::Mask< D, RField<D>, FieldIo<D> >
   {};

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      extern template class Mask< 1, RField<1>, Rpg::FieldIo<1> >;
      extern template class Mask< 2, RField<2>, Rpg::FieldIo<2> >;
      extern template class Mask< 3, RField<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      extern template class Mask<1>;
      extern template class Mask<2>;
      extern template class Mask<3>;
   }
} 
#endif
