#ifndef RPC_C_FIELDS_H
#define RPC_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>       // base class template
#include <rpc/system/Types.h>       // template parameter

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class FieldIo;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A container for c fields stored in both basis and r-grid format.
   *
   * Specializations of this template with D =1, 2, and 3 are derived
   * from specializations of the base class template Rp::CFields, and
   * inherit their public interface and all of their source code from
   * this base class. 
   *
   * \see Rp::CFields
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class CFields : public Rp::CFields<D, Types<D> >
   {};

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      extern template class CFields<1, Rpc::Types<1> >;
      extern template class CFields<2, Rpc::Types<2> >;
      extern template class CFields<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class CFields<1>;
      extern template class CFields<2>;
      extern template class CFields<3>;
   }
}
#endif
