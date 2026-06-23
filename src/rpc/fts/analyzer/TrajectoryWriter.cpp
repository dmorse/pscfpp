/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"                    // header
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/WFields.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <rp/fts/analyzer/TrajectoryWriter.tpp>  // implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   TrajectoryWriter<D>::TrajectoryWriter(Simulator<D>& simulator, 
                                         Rp::System<D, Rpc::Types<D> >& system)
    : Rp::TrajectoryWriter< D, Types<D> > (simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryWriter<1, Rpc::Types<1> >;
      template class TrajectoryWriter<2, Rpc::Types<2> >;
      template class TrajectoryWriter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class TrajectoryWriter<1>;
      template class TrajectoryWriter<2>;
      template class TrajectoryWriter<3>;
   }
}
