/*
* PSCF - MixtureModifier Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureModifier.h"
#include <prdc/solvers/MixtureModifierPrdc.tpp>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/solvers/Block.h>

namespace Pscf {
   namespace Prdc { 
      template class MixtureModifierPrdc< Rpg::Mixture<1> >;
      template class MixtureModifierPrdc< Rpg::Mixture<2> >;
      template class MixtureModifierPrdc< Rpg::Mixture<3> >;
   }
   namespace Rpg { 
      template class MixtureModifier<1>;
      template class MixtureModifier<2>;
      template class MixtureModifier<3>;
   }
}
