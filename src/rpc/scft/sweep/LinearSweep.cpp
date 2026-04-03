/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearSweep.h"
#include <rpc/system/System.h>

#include <rp/scft/sweep/LinearSweep.tpp>

namespace Pscf {
namespace Rpc {

   // Constructor
   template <int D>
   LinearSweep<D>::LinearSweep(System<D>& system)
    : Rp::LinearSweep<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LinearSweep< 1, Rpc::Types<1> >;
      template class LinearSweep< 2, Rpc::Types<2> >;
      template class LinearSweep< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class LinearSweep<1>;
      template class LinearSweep<2>;
      template class LinearSweep<3>;
   }
}
