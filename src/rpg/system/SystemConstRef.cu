/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRef.h"
#include <prdc/system/SystemConstRefTmpl.tpp>

namespace Pscf {
   namespace Prdc {
      template class SystemConstRefTmpl< Rpg::System<1> >;
      template class SystemConstRefTmpl< Rpg::System<2> >;
      template class SystemConstRefTmpl< Rpg::System<3> >;
   }
   namespace Rpg {
      template class SystemConstRef<1>;
      template class SystemConstRef<2>;
      template class SystemConstRef<3>;
   } 
} 
