/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"                    // header
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>

#include <rp/fts/analyzer/ChiDerivative.tpp>  // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   ChiDerivative<D>::ChiDerivative(Simulator<D>& simulator,
                                   System<D>& system)
    : Rp::ChiDerivative< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ChiDerivative< 1, Rpg::Types<1> >;
      template class ChiDerivative< 2, Rpg::Types<2> >;
      template class ChiDerivative< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class ChiDerivative<1>;
      template class ChiDerivative<2>;
      template class ChiDerivative<3>;
   }
}
