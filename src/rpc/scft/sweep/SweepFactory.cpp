/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"                   // header
#include <rpc/scft/sweep/LinearSweep.h>
#include <rp/scft/sweep/SweepFactory.tpp>   // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepFactory<1, Rpc::Types<1> >;
      template class SweepFactory<2, Rpc::Types<2> >;
      template class SweepFactory<3, Rpc::Types<3> >;
   }
}
