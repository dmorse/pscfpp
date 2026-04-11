/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRef.h"
#include <rp/system/SystemConstRef.tpp>

// Explicit instantiations definitions
namespace Pscf {
   namespace Rp {
      template class SystemConstRef< Rpg::System<1> >;
      template class SystemConstRef< Rpg::System<2> >;
      template class SystemConstRef< Rpg::System<3> >;
   }
   namespace Rpg {
      template class SystemConstRef<1>;
      template class SystemConstRef<2>;
      template class SystemConstRef<3>;
   } 
} 
