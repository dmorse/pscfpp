#ifndef RPG_BD_MOVE_H
#define RPG_BD_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/BdMove.h>     // base class template
#include <pscf/cuda/Cuda.h>             // base class argument 
#include <rpg/fts/montecarlo/McMove.h>    // indirect base class
#include <prdc/field/cuda/RField.h>       // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdMove<1, CudaTp<1> >;
      extern template class BdMove<2, CudaTp<2> >;
      extern template class BdMove<3, CudaTp<3> >;
   }
}
#endif
