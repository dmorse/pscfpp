/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"                   // header
#include <rp/scft/sweep/SweepFactory.tpp>   // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepFactory<1,CPT>;
      template class SweepFactory<2,CPT>;
      template class SweepFactory<3,CPT>;
   }
}
