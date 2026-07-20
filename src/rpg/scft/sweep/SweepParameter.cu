/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepParameter.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/MixtureModifier.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/solvers/Block.h>
#include <rpg/field/Domain.h>

#include <rp/scft/sweep/SweepParameter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepParameter<1, CudaTp<1> >;
      template class SweepParameter<2, CudaTp<2> >;
      template class SweepParameter<3, CudaTp<3> >;
      
   }
}
