/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"                    // header

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <rpg/field/WFields.h>

#include <rp/fts/analyzer/TrajectoryWriter.tpp>  // implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   TrajectoryWriter<D>::TrajectoryWriter(Rp::Simulator<D, Rpg::Types<D> >& simulator, 
                                         Rp::System<D, Rpg::Types<D> >& system)
    : Rp::TrajectoryWriter< D, Types<D> > (simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryWriter<1, Rpg::Types<1> >;
      template class TrajectoryWriter<2, Rpg::Types<2> >;
      template class TrajectoryWriter<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class TrajectoryWriter<1>;
      template class TrajectoryWriter<2>;
      template class TrajectoryWriter<3>;
   }
}
