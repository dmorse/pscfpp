/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Solvent.tpp>
#include <pscf/backend/cpp/CPT.h>

namespace Pscf {
   namespace Rp {
      template class Solvent<1,CPT>;
      template class Solvent<2,CPT>;
      template class Solvent<3,CPT>;
   }
}
