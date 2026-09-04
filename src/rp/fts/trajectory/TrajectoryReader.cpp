/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"
#include <pscf/backend/CPT.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class TrajectoryReader<1,CPT>;
      template class TrajectoryReader<2,CPT>;
      template class TrajectoryReader<3,CPT>;
   }
}
