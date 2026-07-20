#ifndef RPC_SYSTEM_CONST_REF_H
#define RPC_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>   // base class template
#include <pscf/cpu/CppTp.h>           // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef<1, CppTp<1> >;
      extern template class SystemConstRef<2, CppTp<2> >;
      extern template class SystemConstRef<3, CppTp<3> >;
   }
}
#endif
