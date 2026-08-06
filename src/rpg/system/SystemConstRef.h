#ifndef RPG_SYSTEM_CONST_REF_H
#define RPG_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>          // template parameter
#include <rp/system/SystemConstRef.h>   // class template

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SystemConstRef<1,CUT>;
      extern template class SystemConstRef<2,CUT>;
      extern template class SystemConstRef<3,CUT>;
   }
}
#endif
