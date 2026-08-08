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
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>

#include <rp/fts/analyzer/TrajectoryWriter.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryWriter<1,CUT>;
      template class TrajectoryWriter<2,CUT>;
      template class TrajectoryWriter<3,CUT>;
   }
}
