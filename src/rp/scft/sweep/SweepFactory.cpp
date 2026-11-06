/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>          // backend identifier class
#include <rp/scft/sweep/SweepFactory.tpp>  // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepFactory<1,CPT>;
      template class SweepFactory<2,CPT>;
      template class SweepFactory<3,CPT>;
   }
}
