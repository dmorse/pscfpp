/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepParameter.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/MixtureModifier.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/solvers/Block.h>
#include <rp/field/Domain.h>

#include <rp/scft/sweep/SweepParameter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepParameter<1,CPT>;
      template class SweepParameter<2,CPT>;
      template class SweepParameter<3,CPT>;
      
   }
}
