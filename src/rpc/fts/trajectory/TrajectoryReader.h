#ifndef RPC_TRAJECTORY_READER_H
#define RPC_TRAJECTORY_READER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/trajectory/TrajectoryReader.h>
#include <pscf/cpu/CppTp.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryReader<1, CppTp<1> >;
      extern template class TrajectoryReader<2, CppTp<2> >;
      extern template class TrajectoryReader<3, CppTp<3> >;
   }
}
#endif
