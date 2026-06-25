/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Ramp.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rp/fts/ramp/Ramp.tpp>     // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Ramp<1, Rpc::Types<1> >;
      template class Ramp<2, Rpc::Types<2> >;
      template class Ramp<3, Rpc::Types<3> >;
   }
}
