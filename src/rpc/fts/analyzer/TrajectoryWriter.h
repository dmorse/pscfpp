#ifndef RPC_TRAJECTORY_WRITER_H
#define RPC_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/TrajectoryWriter.h> // base class template
#include <pscf/cpu/Cpp.h>                 // base class argument
#include <rpc/fts/analyzer/Analyzer.h>        // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryWriter<1, Cpp<1> >;
      extern template class TrajectoryWriter<2, Cpp<2> >;
      extern template class TrajectoryWriter<3, Cpp<3> >;
   }
}
#endif
