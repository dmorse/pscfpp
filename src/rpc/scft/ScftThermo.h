#ifndef RPC_SCFT_THERMO_H
#define RPC_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.h>          // base class template
#include <rpc/system/Types.h>            // template argument
#include <rpc/system/SystemConstRef.h>   // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ScftThermo<1, Rpc::Types<1> >;
      extern template class Rp::ScftThermo<2, Rpc::Types<2> >;
      extern template class Rp::ScftThermo<3, Rpc::Types<3> >;
   }
} 
#endif
