/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimState.h"
#include <prdc/field/cuda/RField.h>
#include <rp/fts/simulator/SimState.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SimState<1, Rpg::Types<1> >;
      template class SimState<2, Rpg::Types<2> >;
      template class SimState<3, Rpg::Types<3> >;
   }
}
