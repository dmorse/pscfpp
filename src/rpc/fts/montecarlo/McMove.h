#ifndef RPC_MC_MOVE_H
#define RPC_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>
#include <rpc/system/Types.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMove<1, Rpc::Types<1> >;
      extern template class McMove<2, Rpc::Types<2> >;
      extern template class McMove<3, Rpc::Types<3> >;
   }
}
#endif
