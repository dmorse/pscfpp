#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <rp/system/System.h>      // template
#include <rpc/system/Types.h>      // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class System<1, Rpc::Types<1> >;
      extern template class System<2, Rpc::Types<1> >;
      extern template class System<3, Rpc::Types<1> >;
   }
}
#endif
