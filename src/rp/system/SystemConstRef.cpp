/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SystemConstRef<1,CPT>;
      template class SystemConstRef<2,CPT>;
      template class SystemConstRef<3,CPT>;
   }
}
