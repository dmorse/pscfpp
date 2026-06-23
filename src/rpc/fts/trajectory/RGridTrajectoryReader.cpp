/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RGridTrajectoryReader.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
//#include <util/misc/ioUtil.h>
//#include <util/misc/FileMaster.h>

#include <rp/fts/trajectory/RGridTrajectoryReader.tpp>

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   RGridTrajectoryReader<D>::RGridTrajectoryReader(Rp::System<D, Rpc::Types<D> >& system)
    : Rp::RGridTrajectoryReader<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class RGridTrajectoryReader<1, Rpc::Types<1> >;
      template class RGridTrajectoryReader<2, Rpc::Types<2> >;
      template class RGridTrajectoryReader<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class RGridTrajectoryReader<1>;
      template class RGridTrajectoryReader<2>;
      template class RGridTrajectoryReader<3>;
   }
}
