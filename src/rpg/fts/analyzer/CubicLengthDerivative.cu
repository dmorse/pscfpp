/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CubicLengthDerivative.h"

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <pscf/cuda/Reduce.h>

#include <rp/fts/analyzer/CubicLengthDerivative.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Constructor.
   template <int D>
   CubicLengthDerivative<D>::CubicLengthDerivative(
                                               Simulator<D>& simulator,
                                               System<D>& system)
    : Rp::CubicLengthDerivative<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class CubicLengthDerivative<1, Rpg::Types<1> >;
      template class CubicLengthDerivative<2, Rpg::Types<2> >;
      template class CubicLengthDerivative<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class CubicLengthDerivative<1>;
      template class CubicLengthDerivative<2>;
      template class CubicLengthDerivative<3>;
   }
}
