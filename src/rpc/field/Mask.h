#ifndef RPC_MASK_H
#define RPC_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Mask.h>       // base class template
#include <rpc/system/Types.h>    // base class template argument
#include <prdc/cpu/RField.h>     // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A field to which the total density is constrained.
   *
   * Please refer to the documentation of the base class template Rp::Mask
   * for more complete API documentation for this class template. The
   * public interface of Rpc::Mask is identical to that of the base class
   * template Rp::Mask.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Mask : public Rp::Mask<D, Types<D> >
   {};

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Mask< 1, Rpc::Types<1> >;
      extern template class Mask< 2, Rpc::Types<2> >;
      extern template class Mask< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Mask<1>;
      extern template class Mask<2>;
      extern template class Mask<3>;
   }
}
#endif
