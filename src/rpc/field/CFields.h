#ifndef RPC_C_FIELDS_H
#define RPC_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>       // base class template
#include <prdc/cpu/RField.h>        // base class member

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
   class CFields : public Rp::CFields<D, RField<D>, FieldIo<D> >
   {};

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      extern template class CFields<1, Cpu::RField<1>, Rpc::FieldIo<1> >;
      extern template class CFields<2, Cpu::RField<2>, Rpc::FieldIo<2> >;
      extern template class CFields<3, Cpu::RField<3>, Rpc::FieldIo<3> >;
   }
   namespace Rpc {
      extern template class CFields<1>;
      extern template class CFields<2>;
      extern template class CFields<3>;
   }
}
#endif
