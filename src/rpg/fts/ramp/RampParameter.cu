/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RampParameter.h"
#include <rpg/fts/simulator/Simulator.h>
#include <rp/system/System.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/MixtureModifier.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/solvers/Block.h>
#include <rp/field/Domain.h>

#include <rp/fts/ramp/RampParameter.tpp>    // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RampParameter<1,CUT>;
      template class RampParameter<2,CUT>;
      template class RampParameter<3,CUT>;
   }
}
