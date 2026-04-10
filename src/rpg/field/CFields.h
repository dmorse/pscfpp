#ifndef RPG_C_FIELD_CONTAINER_H
#define RPG_C_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>     // base class template
#include <rpg/field/FieldIo.h>    // base class template argument
#include <prdc/cuda/RField.h>     // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * A list of c fields stored in both basis and r-grid format.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::CFields, 
   * and inherit their public interface and entire implementation from 
   * this base class.  
   *
   * \see Rp::CFields
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class CFields : public Rp::CFields<D, RField<D>, FieldIo<D> >
   {};

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      extern template class Rp::CFields<1, RField<1>, Rpg::FieldIo<1> >;
      extern template class Rp::CFields<2, RField<2>, Rpg::FieldIo<2> >;
      extern template class Rp::CFields<3, RField<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      extern template class CFields<1>;
      extern template class CFields<2>;
      extern template class CFields<3>;
   }
}
#endif
