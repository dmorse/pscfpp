#ifndef RPG_TRAJECTORY_READER_H
#define RPG_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/TrajectoryReader.h>
#include <rpg/system/Types.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryReader<1, Rpg::Types<1> >;
      extern template class TrajectoryReader<2, Rpg::Types<2> >;
      extern template class TrajectoryReader<3, Rpg::Types<3> >;
   }
}
#endif
