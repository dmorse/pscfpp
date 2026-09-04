/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
#include "RGridTrajectoryReader.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#endif

#include <pscf/backend/CUT.h>
#include <rp/fts/trajectory/RGridTrajectoryReader.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class RGridTrajectoryReader<1,CUT>;
      template class RGridTrajectoryReader<2,CUT>;
      template class RGridTrajectoryReader<3,CUT>;
   }
}
