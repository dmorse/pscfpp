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
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#endif

#include <pscf/backends/CPT.h>
#include <rp/fts/trajectory/RGridTrajectoryReader.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class RGridTrajectoryReader<1,CPT>;
      template class RGridTrajectoryReader<2,CPT>;
      template class RGridTrajectoryReader<3,CPT>;
   }
}
