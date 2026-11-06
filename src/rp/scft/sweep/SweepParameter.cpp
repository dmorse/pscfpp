/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <rp/scft/sweep/SweepParameter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepParameter<1,CPT>;
      template class SweepParameter<2,CPT>;
      template class SweepParameter<3,CPT>;
      
   }
}
