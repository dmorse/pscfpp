#ifndef RPC_C_FIELDS_H
#define RPC_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>       // base class template
#include <rpc/system/Types.h>       // template parameter

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CFields<1, Rpc::Types<1> >;
      extern template class CFields<2, Rpc::Types<2> >;
      extern template class CFields<3, Rpc::Types<3> >;
   }
}
#endif
