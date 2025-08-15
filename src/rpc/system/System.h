#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <prdc/system/SystemReal.h>      // base class template
#include <rpc/system/Types.h>            // base class template param
#include <rpc/solvers/Mixture.h>         // member
#include <rpc/field/Domain.h>            // member
#include <rpc/field/WFieldContainer.h>   // member
#include <rpc/field/CFieldContainer.h>   // member
#include <rpc/field/Mask.h>              // member

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /**
   * Main class, representing a complete physical system.
   *
   * This class is essentially a named partial specialization of the base
   * class template SystemReal, and has the same public interface as this 
   * base class.  See the documentation of the base class for the class
   * interface.
   *
   * \ingroup Pscf_Rpc_Module
   */
   template <int D>
   class System : public SystemReal< D, Types<D> >
   {
   public:

      /**
      * Constructor.
      */
      System();

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
   extern template class SystemReal<1, Rpc::Types<1> >;
   extern template class SystemReal<2, Rpc::Types<1> >;
   extern template class SystemReal<3, Rpc::Types<1> >;
   #endif

} // namespace Prdc 
} // namespace Pscf
#endif
