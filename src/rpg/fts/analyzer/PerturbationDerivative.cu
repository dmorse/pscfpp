/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"

#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <rp/fts/analyzer/PerturbationDerivative.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationDerivative<1,CUT>;
      template class PerturbationDerivative<2,CUT>;
      template class PerturbationDerivative<3,CUT>;
   }
}
