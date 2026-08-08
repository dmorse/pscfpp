/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <rp/fts/analyzer/ConcentrationDerivative.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationDerivative<1,CUT>;
      template class ConcentrationDerivative<2,CUT>;
      template class ConcentrationDerivative<3,CUT>;
   }
}
