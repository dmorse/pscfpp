/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRef.h"
#include <prdc/system/SystemConstRefReal.tpp>

namespace Pscf {
   namespace Prdc {
      template class SystemConstRefReal< Rpc::System<1> >;
      template class SystemConstRefReal< Rpc::System<2> >;
      template class SystemConstRefReal< Rpc::System<3> >;
   }
   namespace Rpc {
      template class SystemConstRef<1>;
      template class SystemConstRef<2>;
      template class SystemConstRef<3>;
   } 
} 
