/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <rp/fts/ramp/RampFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RampFactory<1,CPT>;
      template class RampFactory<2,CPT>;
      template class RampFactory<3,CPT>;
   }
}
