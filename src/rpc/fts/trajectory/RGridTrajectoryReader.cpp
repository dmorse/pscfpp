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

#include <rp/fts/trajectory/RGridTrajectoryReader.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class RGridTrajectoryReader<1,CPT>;
      template class RGridTrajectoryReader<2,CPT>;
      template class RGridTrajectoryReader<3,CPT>;
   }
}
