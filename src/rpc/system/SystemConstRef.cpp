/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRef.h"
#include <rpc/system/System.h>
#include <rpc/field/Domain.h>

#include <rp/system/SystemConstRef.tpp>

// Explicit initialization definitions
namespace Pscf {
   namespace Rp {
      template class SystemConstRef<1, Rpc::Types<1> >;
      template class SystemConstRef<2, Rpc::Types<2> >;
      template class SystemConstRef<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class SystemConstRef<1>;
      template class SystemConstRef<2>;
      template class SystemConstRef<3>;
   } 
} 
