#ifndef RPC_DOMAIN_H
#define RPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Domain.h>     // base class template
#include <pscf/cpu/Cpp.h>    // base class template argument

// Explicit instantiation declarations 
namespace Pscf {
   namespace Rp {
      extern template class Domain<1, Cpp<1> >;
      extern template class Domain<2, Cpp<2> >;
      extern template class Domain<3, Cpp<3> >;
   } 
} 
#endif
