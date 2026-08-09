#ifndef RPC_RGRID_TRAJECTORY_READER_H
#define RPC_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/RGridTrajectoryReader.h> // direct base class 
#include <pscf/backends/CPT.h>                        // base class argument
#include <prdc/field/cpu/RField.h>                   // base class member
#include <rp/fts/trajectory/TrajectoryReader.h>     // indirect base class 

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RGridTrajectoryReader<1,CPT>;
      extern template class RGridTrajectoryReader<2,CPT>;
      extern template class RGridTrajectoryReader<3,CPT>;
   }
}
#endif
