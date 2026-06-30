/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"

#if 0
namespace Pscf {
namespace Rpc {

   /*
   * Constructor.
   */
   template <int D>
   Rp::TrajectoryReader<D, Rpc::Types<D> >::TrajectoryReader(Rp::System<D, Rpc::Types<D> >& system)
    : Rp::TrajectoryReader<D, Types<D> >(system)
   {}

}
}
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class TrajectoryReader<1, Rpc::Types<1> >;
      template class TrajectoryReader<2, Rpc::Types<2> >;
      template class TrajectoryReader<3, Rpc::Types<3> >;
   }
}
