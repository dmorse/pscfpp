#ifndef RPC_RGRID_TRAJECTORY_READER_H
#define RPC_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/RGridTrajectoryReader.h> // direct base class 
#include <pscf/cpu/CppTp.h>                        // base class argument
#include <prdc/field/cpu/RField.h>                   // base class member
#include <rpc/fts/trajectory/TrajectoryReader.h>     // indirect base class 

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RGridTrajectoryReader<1, CppTp<1> >;
      extern template class RGridTrajectoryReader<2, CppTp<2> >;
      extern template class RGridTrajectoryReader<3, CppTp<3> >;
   }
}
#endif
