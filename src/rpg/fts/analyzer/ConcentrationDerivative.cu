/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"                    // header
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <rp/fts/analyzer/ConcentrationDerivative.tpp>  // implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   ConcentrationDerivative<D>::ConcentrationDerivative(
                                   Simulator<D>& simulator,
                                   Rp::System<D, Rpg::Types<D> >& system)
    : Rp::ConcentrationDerivative< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationDerivative< 1, Rpg::Types<1> >;
      template class ConcentrationDerivative< 2, Rpg::Types<2> >;
      template class ConcentrationDerivative< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class ConcentrationDerivative<1>;
      template class ConcentrationDerivative<2>;
      template class ConcentrationDerivative<3>;
   }
}
