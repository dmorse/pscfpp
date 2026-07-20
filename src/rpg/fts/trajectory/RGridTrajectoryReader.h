#ifndef RPG_RGRID_TRAJECTORY_READER_H
#define RPG_RGRID_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/RGridTrajectoryReader.h> // direct base class 
#include <pscf/cuda/Cuda.h>                        // base class argument
#include <prdc/field/cuda/RField.h>                  // base class member
#include <rpg/fts/trajectory/TrajectoryReader.h>     // indirect base class 

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RGridTrajectoryReader<1, CudaTp<1> >;
      extern template class RGridTrajectoryReader<2, CudaTp<2> >;
      extern template class RGridTrajectoryReader<3, CudaTp<3> >;
   }
}
#endif
