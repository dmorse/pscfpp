/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RGridTrajectoryReader.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <rp/fts/trajectory/RGridTrajectoryReader.tpp>

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   RGridTrajectoryReader<D>::RGridTrajectoryReader(System<D>& system)
    : Rp::RGridTrajectoryReader<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class RGridTrajectoryReader<1, Rpg::Types<1> >;
      template class RGridTrajectoryReader<2, Rpg::Types<2> >;
      template class RGridTrajectoryReader<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class RGridTrajectoryReader<1>;
      template class RGridTrajectoryReader<2>;
      template class RGridTrajectoryReader<3>;
   }
}
