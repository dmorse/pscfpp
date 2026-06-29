#ifndef RPC_SOLVENT_H
#define RPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Solvent.h>    // base class template
#include <rpc/system/Types.h>      // base class template parameter
#include <prdc/field/cpu/RField.h>       // member of base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Solvent<1, Rpc::Types<1> >;
      extern template class Solvent<2, Rpc::Types<2> >;
      extern template class Solvent<3, Rpc::Types<3> >;
   }
}
#endif
