/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearRamp.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rp/fts/ramp/LinearRamp.tpp>     // base class implementation


// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LinearRamp<1, Rpc::Types<1> >;
      template class LinearRamp<2, Rpc::Types<2> >;
      template class LinearRamp<3, Rpc::Types<3> >;
   }
}
