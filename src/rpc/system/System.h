#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <rp/system/System.h>      // base class template
#include <pscf/cpu/CppTp.h>      // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class System<1, CppTp<1> >;
      extern template class System<2, CppTp<1> >;
      extern template class System<3, CppTp<1> >;
   }
}
#endif
