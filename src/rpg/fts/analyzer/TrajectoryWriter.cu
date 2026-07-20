/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryWriter.h"

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <rpg/field/WFields.h>

#include <rp/fts/analyzer/TrajectoryWriter.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryWriter<1, CudaTp<1> >;
      template class TrajectoryWriter<2, CudaTp<2> >;
      template class TrajectoryWriter<3, CudaTp<3> >;
   }
}
