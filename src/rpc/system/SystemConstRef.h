#ifndef RPC_SYSTEM_CONST_REF_H
#define RPC_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>   // base class template
#include <rpc/system/System.h>          // template parameter

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Const access to a System<D>.
   *
   * Specializations of this template with D=1, 2, and 3 are derived
   * from corresponding specializations of base class template
   * Rp::SystemConstRef, and inherit their public interface and almost
   * all of their source code from this base class.
   *
   * \see Rp::SystemConstRef
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class SystemConstRef : public Rp::SystemConstRef< System<D> >
   {
   public:

      /**
      * Default constructor.
      */
      SystemConstRef()
       : Rp::SystemConstRef< System<D> > ()
      {};

      /**
      * Constructor.
      *
      * \param system  System<D> object to which this refers.
      */
      SystemConstRef(System<D> const & system)
       : Rp::SystemConstRef< System<D> > (system)
      {};

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef< Rpc::System<1> >;
      extern template class SystemConstRef< Rpc::System<2> >;
      extern template class SystemConstRef< Rpc::System<3> >;
   }
   namespace Rpc {
      extern template class SystemConstRef<1>;
      extern template class SystemConstRef<2>;
      extern template class SystemConstRef<3>;
   }
}
#endif
