#ifndef RPC_BASIS_FIELD_STATE_H
#define RPC_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/BasisFieldState.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::BasisFieldState, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::BasisFieldState
   * \ingroup Rpc_Scft_Sweep_Module
   */
   template <int D>
   class BasisFieldState
    : public Rp::BasisFieldState<D, Types<D> >
   {
   public:

      /**
      * Default constructor.
      */
      BasisFieldState();

      /**
      * Constructor, create association with a parent system.
      *
      * \param system associated parent system
      */
      BasisFieldState(Rp::System<D, Rpc::Types<D> >& system);

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BasisFieldState< 1, Rpc::Types<1> >;
      extern template class BasisFieldState< 2, Rpc::Types<2> >;
      extern template class BasisFieldState< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class BasisFieldState<1>;
      extern template class BasisFieldState<2>;
      extern template class BasisFieldState<3>;
   }
}
#endif
