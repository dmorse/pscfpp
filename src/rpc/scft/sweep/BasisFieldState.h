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

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BasisFieldState< 1, Rpc::Types<1> >;
      extern template class BasisFieldState< 2, Rpc::Types<2> >;
      extern template class BasisFieldState< 3, Rpc::Types<3> >;
   }
}
#endif
