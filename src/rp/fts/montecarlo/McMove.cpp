/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <rp/fts/montecarlo/McMove.tpp>     // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMove<1,CPT>;
      template class McMove<2,CPT>;
      template class McMove<3,CPT>;
   }
}
