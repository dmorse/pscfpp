/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReaderFactory.h"

// Subclasses of TrajectoryReader
#include <rpg/fts/trajectory/RGridTrajectoryReader.h>

#include <rp/fts/trajectory/TrajectoryReaderFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryReaderFactory<1, CudaTp<1> >;
      template class TrajectoryReaderFactory<2, CudaTp<2> >;
      template class TrajectoryReaderFactory<3, CudaTp<3> >;
   }
}
