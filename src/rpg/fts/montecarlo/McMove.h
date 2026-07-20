#ifndef RPG_MC_MOVE_H
#define RPG_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>
#include <pscf/cuda/Cuda.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMove<1, CudaTp<1> >;
      extern template class McMove<2, CudaTp<2> >;
      extern template class McMove<3, CudaTp<3> >;
   }
}
#endif
