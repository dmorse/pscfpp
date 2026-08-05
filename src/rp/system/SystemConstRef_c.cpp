/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.tpp>   // class implementation
#include <pscf/backends/CPT.h>            // backend type class

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SystemConstRef<1,CPT>;
      template class SystemConstRef<2,CPT>;
      template class SystemConstRef<3,CPT>;
   }
}
