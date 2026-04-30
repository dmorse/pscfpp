#ifndef RPC_SYSTEM_CONST_REF_H
#define RPC_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>   // base class template
#include <rpc/system/Types.h>           // base class parameter

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Const access to a System.
   *
   * Specializations of this template with D=1, 2, and 3 are derived
   * from corresponding specializations of the base class template
   * Rp::SystemConstRef, and inherit their public interface and almost
   * all of their source code from this base class.
   *
   * \see Rp::SystemConstRef
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class SystemConstRef : public Rp::SystemConstRef<D, Types<D> >
   {
   public:

      /**
      * Default constructor.
      */
      SystemConstRef()
       : Rp::SystemConstRef< D, Types<D> > ()
      {};

      /**
      * Constructor.
      *
      * \param system  System<D> object to which this refers.
      */
      SystemConstRef(System<D> const & system)
       : Rp::SystemConstRef<D, Types<D> > (system)
      {};

      /**
      * Destructor.
      */
      virtual ~SystemConstRef() = default;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef<1, Rpc::Types<1> >;
      extern template class SystemConstRef<2, Rpc::Types<2> >;
      extern template class SystemConstRef<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class SystemConstRef<1>;
      extern template class SystemConstRef<2>;
      extern template class SystemConstRef<3>;
   }
}
#endif
