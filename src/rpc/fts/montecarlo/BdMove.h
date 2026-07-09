#ifndef RPC_BD_MOVE_H
#define RPC_BD_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/BdMove.h>     // base class template
#include <rpc/system/Types.h>             // base class argument 
#include <rpc/fts/montecarlo/McMove.h>    // indirect base class
#include <prdc/field/cpu/RField.h>        // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdMove<1, Rpc::Types<1> >;
      extern template class BdMove<2, Rpc::Types<2> >;
      extern template class BdMove<3, Rpc::Types<3> >;
   }
}
#endif
