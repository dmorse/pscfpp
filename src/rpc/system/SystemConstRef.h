#ifndef RPC_SYSTEM_CONST_REF_H
#define RPC_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/system/SystemConstRefTmpl.h>   // base class template
#include <rpc/system/System.h>                // template parameter

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Const access to a System<D>.
   *
   * See documentation of the base class template Prdc::SystemConstRefTmpl
   * for documentation of the public interface.
   *
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class SystemConstRef : public SystemConstRefTmpl< System<D> >
   {
   public:

      /// Alias for base class
      using Base = SystemConstRefTmpl< System<D> >;

      /**
      * Default constructor.
      */
      SystemConstRef()
       : Base()
      {};

      /**
      * Constructor.
      * 
      * \param system  System<D> object to which this refers.
      */
      SystemConstRef(System<D> const & system)
       : Base(system)
      {};

   };

   // Suppress implicit instantiation
   extern template class SystemConstRef<1>;
   extern template class SystemConstRef<2>;
   extern template class SystemConstRef<3>;

} // namespace Rpc

namespace Prdc {
   // Suppress implicit instantiation of base class
   extern template class SystemConstRefTmpl< Rpc::System<1> >;
   extern template class SystemConstRefTmpl< Rpc::System<2> >;
   extern template class SystemConstRefTmpl< Rpc::System<3> >;
} // namespace Rpc

} // namespace Pscf
#endif
