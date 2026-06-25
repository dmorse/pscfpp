/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"
#include <rpc/fts/simulator/SimState.h>
#include <prdc/cpu/RField.h>

#include <rp/fts/montecarlo/McMove.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMove<1, Rpc::Types<1> >;
      template class McMove<2, Rpc::Types<2> >;
      template class McMove<3, Rpc::Types<3> >;
   }
}
