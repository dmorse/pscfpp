/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>
#include <rp/field/CFields.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class CFields<1,CUT>;
      template class CFields<2,CUT>;
      template class CFields<3,CUT>;
   }
}
