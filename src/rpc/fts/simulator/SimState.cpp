/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimState.h"
#include <prdc/cpu/RField.h>
#include <rp/fts/simulator/SimState.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template struct SimState<1, Rpc::Types<1> >;
      template struct SimState<2, Rpc::Types<2> >;
      template struct SimState<3, Rpc::Types<3> >;
   }
}
