/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationDerivative.h"                    // header
#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <rp/fts/analyzer/ConcentrationDerivative.tpp>  // implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   ConcentrationDerivative<D>::ConcentrationDerivative(
                                   Simulator<D>& simulator,
                                   System<D>& system)
    : Rp::ConcentrationDerivative< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationDerivative< 1, Rpc::Types<1> >;
      template class ConcentrationDerivative< 2, Rpc::Types<2> >;
      template class ConcentrationDerivative< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class ConcentrationDerivative<1>;
      template class ConcentrationDerivative<2>;
      template class ConcentrationDerivative<3>;
   }
}
