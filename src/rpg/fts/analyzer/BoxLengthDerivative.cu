/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BoxLengthDerivative.h"

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <pscf/cuda/Reduce.h>

#include <rp/fts/analyzer/BoxLengthDerivative.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BoxLengthDerivative<D>::BoxLengthDerivative(Simulator<D>& simulator,
                                               System<D>& system)
    : Rp::BoxLengthDerivative<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class BoxLengthDerivative<1, Rpg::Types<1> >;
      template class BoxLengthDerivative<2, Rpg::Types<2> >;
      template class BoxLengthDerivative<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class BoxLengthDerivative<1>;
      template class BoxLengthDerivative<2>;
      template class BoxLengthDerivative<3>;
   }
}
