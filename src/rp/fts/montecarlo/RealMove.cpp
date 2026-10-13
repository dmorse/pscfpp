
/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/RealMove.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RealMove<1,CPT>;
      template class RealMove<2,CPT>;
      template class RealMove<3,CPT>;
   }
}
