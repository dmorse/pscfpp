/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/BdMove.tpp> 
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdMove<1,CPT>;
      template class BdMove<2,CPT>;
      template class BdMove<3,CPT>;
   }
}
