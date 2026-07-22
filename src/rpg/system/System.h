#ifndef RPG_SYSTEM_H
#define RPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <rp/system/System.h>       // base class template
#include <pscf/backends/CUT.h>       // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class System<1,CUT>;
      extern template class System<2,CUT>;
      extern template class System<3,CUT>;
   }
}
#endif
