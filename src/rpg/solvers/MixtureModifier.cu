/*
* PSCF - MixtureModifier Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureModifier.h"
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/solvers/Block.h>

#include <rp/solvers/MixtureModifier.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MixtureModifier< 1, CudaTp<1> >;
      template class MixtureModifier< 2, CudaTp<2> >;
      template class MixtureModifier< 3, CudaTp<3> >;
   }
}
