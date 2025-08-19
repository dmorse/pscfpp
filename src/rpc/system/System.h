#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <prdc/system/SystemTmpl.h>      // base class template
#include <rpc/system/Types.h>            // base class template param
#include <rpc/solvers/Mixture.h>         // member
#include <rpc/field/Domain.h>            // member
#include <rpc/field/WFields.h>   // member
#include <rpc/field/CFieldContainer.h>   // member
#include <rpc/field/Mask.h>              // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Main class, representing a complete physical system.
   *
   * This class is derived from a partial specialization of the class
   * template Prdc::SystemTmpl, and has the same public interface as 
   * its base class.  See the documentation of this base class template 
   * for details.
   *
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class System : public SystemTmpl< D, Types<D> >
   {
   public:

      /**
      * Constructor.
      */
      System();

      // Suppress compiler-generated member functions
      System(System<D> const &) = delete;
      System<D>& operator = (System<D> const &) = delete;

   };

   #ifndef RPC_SYSTEM_CPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Rpc
namespace Prdc {

   #ifndef RPC_SYSTEM_CPP
   // Suppress implicit instantiation
   extern template class SystemTmpl<1, Rpc::Types<1> >;
   extern template class SystemTmpl<2, Rpc::Types<1> >;
   extern template class SystemTmpl<3, Rpc::Types<1> >;
   #endif

} // namespace Prdc 
} // namespace Pscf
#endif
