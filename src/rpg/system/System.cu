#ifndef RPG_SYSTEM_CU
#define RPG_SYSTEM_CU

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>    // backend type class
#include <rp/system/System.tpp>   // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class System<1,CUT>;
      template class System<2,CUT>;
      template class System<3,CUT>;
   }
}
#endif
