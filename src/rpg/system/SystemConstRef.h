#ifndef RPG_SYSTEM_CONST_REF_H
#define RPG_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>   // base class template
#include <rpg/system/Types.h>           // base class parameter

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef<1, Rpg::Types<1> >;
      extern template class SystemConstRef<2, Rpg::Types<2> >;
      extern template class SystemConstRef<3, Rpg::Types<3> >;
   }
}
#endif
