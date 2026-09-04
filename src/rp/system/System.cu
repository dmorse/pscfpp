/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/System.tpp>  // class implementation
#include <pscf/backend/CUT.h>   // backend type class

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class System<1,CUT>;
      template class System<2,CUT>;
      template class System<3,CUT>;
   }
}
