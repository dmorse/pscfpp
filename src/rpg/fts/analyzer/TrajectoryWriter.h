#ifndef RPG_TRAJECTORY_WRITER_H
#define RPG_TRAJECTORY_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/TrajectoryWriter.h> // base template
#include <pscf/cuda/CudaTp.h>                 // base argument
#include "Analyzer.h"                         // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class TrajectoryWriter<1, CudaTp<1> >;
      extern template class TrajectoryWriter<2, CudaTp<2> >;
      extern template class TrajectoryWriter<3, CudaTp<3> >;
   }
}
#endif
