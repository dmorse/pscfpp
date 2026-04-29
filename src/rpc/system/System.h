#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <rp/system/System.h>      // base class template
#include <rpc/system/Types.h>      // base class template argument

// Forward declarations
namespace Pscf {
   namespace Rpc {
      template <int D> class WFields;
      template <int D> class CFields;
      template <int D> class Mask;
   }
}
namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * A complete physical system.
   *
   * Specializations of this template with D=1, 2 or 3 are derived from
   * corresponding specializations of the base class template Rp::System,
   * and each have the same public interface as this base class.
   *
   * \see Rp::System
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class System : public Rp::System< D, Types<D> >
   {
   public:

      /**
      * Constructor.
      */
      System();

      // Deleted functions
      System(System<D> const &) = delete;
      System<D>& operator = (System<D> const &) = delete;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class System<1, Rpc::Types<1> >;
      extern template class System<2, Rpc::Types<1> >;
      extern template class System<3, Rpc::Types<1> >;
   }
   namespace Rpc {
      extern template class System<1>;
      extern template class System<2>;
      extern template class System<3>;
   }
}
#endif
