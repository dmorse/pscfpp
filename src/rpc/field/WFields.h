#ifndef RPC_W_FIELDS_H
#define RPC_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/WFields.h>    // base class template
#include <rpc/system/Types.h>    // base class member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * Specializations of this template with D =1, 2, and 3 are derived
   * from specializations of the base class template Rp::WFields, and
   * inherit their public interface and all of their source code
   * from this base class. See the documentation for this base class
   * template for details.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class WFields : public Rp::WFields<D, Types<D> >
   {};

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      extern template class WFields< 1, Rpc::Types<1> >;
      extern template class WFields< 2, Rpc::Types<2> >;
      extern template class WFields< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class WFields<1>;
      extern template class WFields<2>;
      extern template class WFields<3>;
   }
}
#endif
